library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(jsonlite)
library(tidyr)

mutation_list_raw_import <- c("D178N", "E200K", "P102L")

# Return the fallback value when the left-hand value is NULL.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Standardise assay labels and keep only the three expected SNV assays.
normalise_assay <- function(x) {
  # Force a consistent case before matching assay names.
  x <- toupper(as.character(x))

  # Correct the recurrent D178N typo before extracting the recognised assay.
  x <- str_replace_all(x, "D1789N", "D178N")
  str_extract(x, "D178N|E200K|P102L")
}

# Clean raw target names into analysis labels, including WT and assay typo fixes.
clean_ddpcr_target <- function(x) {
  # Remove exported probe/channel suffixes from target labels.
  out <- as.character(x)
  out <- str_replace_all(out, "-mut|_FAM1|_VIC2", "")

  # Correct the recurrent D178N typo and map PRNP to the WT reference label.
  out <- str_replace_all(out, "D1789N", "D178N")
  out[out == "PRNP"] <- "WT"
  out
}

# Interpret metadata values that may be stored as logical TRUE or the string "true".
is_true <- function(x) {
  isTRUE(x) || identical(tolower(as.character(x)), "true")
}

# Read one scalar field from a target metadata object, or return a default.
target_value <- function(target, key, default = NA_character_) {
  # Metadata fields may be absent or empty for some target formats.
  value <- target[[key]]
  if (is.null(value) || length(value) == 0L) {
    return(default)
  }

  # Keep a scalar string even if the metadata field is represented as a list.
  as.character(value[[1]])
}

# Return the target name from either TargetName or, if absent, Name.
target_name <- function(target) {
  # Prefer the explicit TargetName field used by newer metadata exports.
  value <- target_value(target, "TargetName")

  # Fall back to Name for older or alternate metadata structures.
  if (is.na(value)) {
    value <- target_value(target, "Name")
  }
  value
}

# Return the ddPCR fluorescence channel recorded for a target metadata object.
target_channel <- function(target) {
  # Newer metadata can store a single dye object with the channel directly.
  dye <- target[["Dye"]]
  if (!is.null(dye) && !is.null(dye[["Channel"]])) {
    return(as.integer(dye[["Channel"]]))
  }

  # Older metadata can store one or more dye objects in a Dyes list.
  dyes <- target[["Dyes"]]
  if (!is.null(dyes)) {
    # Some metadata stores dyes as a list; keep the first non-missing channel.
    channels <- purrr::map_int(
      dyes,
      # Pull the channel number from one dye entry.
      function(dye_entry) {
        if (is.null(dye_entry) || is.null(dye_entry[["Channel"]])) {
          return(NA_integer_)
        }
        as.integer(dye_entry[["Channel"]])
      }
    )
    channels <- channels[!is.na(channels)]
    if (length(channels) > 0L) {
      return(channels[[1]])
    }
  }

  # Return a missing channel if neither metadata layout contains one.
  NA_integer_
}

# Count droplets assigned to one metadata cluster, treating absent droplet lists as zero.
cluster_droplet_count <- function(cluster) {
  droplets <- cluster[["Droplets"]]
  if (is.null(droplets)) {
    return(0L)
  }
  length(droplets)
}

lambda_cap <- 1000

# Convert reference and mutant lambdas into percentage mutant fractional abundance.
fractional_abundance_from_lambdas <- function(ref_lambda, mut_lambda) {
  # Missing lambdas cannot produce an interpretable fractional abundance.
  if (is.na(ref_lambda) || is.na(mut_lambda)) {
    return(NA_real_)
  }

  # The denominator is the total Poisson-corrected concentration signal.
  denominator <- ref_lambda + mut_lambda
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  100 * mut_lambda / denominator
}

# Estimate lambda directly from a positive partition count and total droplets.
lambda_from_positive_count <- function(positive, total) {
  # Reject impossible counts before applying the Poisson transform.
  if (is.na(positive) || is.na(total) || total <= 0 || positive < 0 || positive > total) {
    return(NA_real_)
  }

  # Zero positives imply zero estimated concentration.
  if (positive <= 0) {
    return(0)
  }

  # Fully positive wells saturate the transform, so cap the estimate.
  if (positive >= total) {
    return(lambda_cap)
  }

  -log1p(-positive / total)
}

# Compute a binomial confidence interval for lambda from positive droplet counts.
binomial_lambda_interval <- function(positive, total, conf.level = 0.95) {
  # Work with numeric scalar counts.
  positive <- as.numeric(positive)
  total <- as.numeric(total)

  # Return an all-missing interval for invalid or impossible inputs.
  if (is.na(positive) || is.na(total) || total <= 0 || positive < 0 || positive > total) {
    return(tibble(lambda = NA_real_, lower = NA_real_, upper = NA_real_))
  }

  # Compute exact binomial bounds on the positive-droplet proportion.
  alpha <- 1 - conf.level
  p_lower <- if (positive <= 0) {
    0
  } else {
    qbeta(alpha / 2, positive, total - positive + 1)
  }
  p_upper <- if (positive >= total) {
    1
  } else {
    qbeta(1 - alpha / 2, positive + 1, total - positive)
  }

  # Transform the point estimate and interval bounds from proportion to lambda.
  tibble(
    lambda = lambda_from_positive_count(positive, total),
    lower = lambda_from_positive_count(total * p_lower, total),
    upper = ifelse(
      p_upper >= 1,
      lambda_cap,
      lambda_from_positive_count(total * p_upper, total)
    )
  )
}

# Compute concentration-ratio confidence limits using Fieller-style interval geometry.
fieller_ratio_interval <- function(
  numerator,
  denominator,
  numerator_lower,
  numerator_upper,
  denominator_lower,
  denominator_upper
) {
  values <- c(
    numerator, denominator,
    numerator_lower, numerator_upper,
    denominator_lower, denominator_upper
  )
  if (any(is.na(values)) || denominator <= 0) {
    return(tibble(ratio_lower = NA_real_, ratio_upper = NA_real_, ratio_unbounded = NA))
  }

  # Find the tangent slopes from the origin to an uncertainty ellipse quadrant.
  tangent_slopes <- function(center_x, center_y, radius_x, radius_y) {
    a <- center_x^2 - radius_x^2
    discriminant <- center_x^2 * radius_y^2 +
      radius_x^2 * center_y^2 -
      radius_x^2 * radius_y^2

    if (a <= 0) {
      return(c(-Inf, Inf))
    }
    if (discriminant < 0) {
      return(c(NA_real_, NA_real_))
    }

    sort(c(
      (center_x * center_y - sqrt(discriminant)) / a,
      (center_x * center_y + sqrt(discriminant)) / a
    ))
  }

  # Calculate lower and upper tangent bounds for the ratio interval.
  lower_slopes <- tangent_slopes(
    center_x = denominator,
    center_y = numerator,
    radius_x = max(denominator_upper - denominator, 0),
    radius_y = max(numerator - numerator_lower, 0)
  )
  upper_slopes <- tangent_slopes(
    center_x = denominator,
    center_y = numerator,
    radius_x = max(denominator - denominator_lower, 0),
    radius_y = max(numerator_upper - numerator, 0)
  )

  # If the ellipse geometry fails, return missing ratio bounds.
  if (any(is.na(c(lower_slopes, upper_slopes)))) {
    return(tibble(ratio_lower = NA_real_, ratio_upper = NA_real_, ratio_unbounded = NA))
  }

  # Keep ratios non-negative and record whether the upper interval is unbounded.
  tibble(
    ratio_lower = max(0, lower_slopes[[1]]),
    ratio_upper = max(0, upper_slopes[[2]]),
    ratio_unbounded = is.infinite(upper_slopes[[2]])
  )
}

# Estimate mutant fractional abundance and confidence limits from ddPCR counts.
# Inputs are paired reference and mutant positive/negative partition counts for
# the same accepted droplet total. The function first converts each target's
# positive-droplet fraction to a Poisson-corrected concentration signal
# (lambda), then estimates the mutant share of the total target signal.
fractional_abundance <- function(
  ref_positive,
  ref_negative,
  mut_positive,
  mut_negative,
  total,
  conf.level = 0.95
) {
  # Coerce all count inputs to numeric scalars before validation. These are
  # usually pooled droplet counts, not raw row-wise strings.
  ref_positive <- as.numeric(ref_positive)
  ref_negative <- as.numeric(ref_negative)
  mut_positive <- as.numeric(mut_positive)
  mut_negative <- as.numeric(mut_negative)
  total <- as.numeric(total)

  # Count pairs must be complete, non-negative, and sum to the same accepted
  # droplet total. Reference and mutant calls are two views of the same well or
  # pooled wells, so each positive+negative pair should exhaust total droplets.
  invalid <- any(is.na(c(ref_positive, ref_negative, mut_positive, mut_negative, total))) ||
    total <= 0 ||
    ref_positive < 0 || mut_positive < 0 ||
    ref_negative < 0 || mut_negative < 0 ||
    ref_positive + ref_negative != total ||
    mut_positive + mut_negative != total

  # Preserve output shape for invalid rows so pmap_dfr() callers can bind
  # results safely and detect the failed estimate from NA fields.
  if (invalid) {
    return(tibble(
      fractional_abundance = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      lambda_ref = NA_real_,
      lambda_mut = NA_real_,
      ratio = NA_real_,
      ratio_low = NA_real_,
      ratio_high = NA_real_,
      ratio_unbounded = NA
    ))
  }

  # Convert each target's positive droplet count to lambda, the
  # Poisson-corrected concentration estimate in molecules per droplet. The
  # interval is computed on the positive-droplet proportion and transformed to
  # the lambda scale inside binomial_lambda_interval().
  ref <- binomial_lambda_interval(ref_positive, total, conf.level = conf.level)
  mut <- binomial_lambda_interval(mut_positive, total, conf.level = conf.level)

  # The denominator for mutant fractional abundance is the total corrected
  # target concentration: reference signal plus mutant signal.
  total_lambda <- ref$lambda + mut$lambda

  # If the concentration ratio is undefined, still return the point estimate if
  # possible. A ratio-based CI needs a positive finite reference lambda.
  if (!is.finite(total_lambda) || total_lambda <= 0 || ref$lambda <= 0) {
    return(tibble(
      fractional_abundance = fractional_abundance_from_lambdas(ref$lambda, mut$lambda),
      ci_low = NA_real_,
      ci_high = NA_real_,
      lambda_ref = ref$lambda,
      lambda_mut = mut$lambda,
      ratio = NA_real_,
      ratio_low = NA_real_,
      ratio_high = NA_real_,
      ratio_unbounded = NA
    ))
  }

  # Dube et al. derive dPCR concentration-ratio confidence intervals for
  # lambda_target / lambda_reference. Here the same ratio geometry is applied
  # to lambda_mut / lambda_wt. Working on the ratio scale is useful because the
  # mutant fraction is a transformed concentration ratio, not an independent
  # binomial proportion of all droplets.
  ratio_ci <- fieller_ratio_interval(
    numerator = mut$lambda,
    denominator = ref$lambda,
    numerator_lower = mut$lower,
    numerator_upper = mut$upper,
    denominator_lower = ref$lower,
    denominator_upper = ref$upper
  )
  ratio <- mut$lambda / ref$lambda

  # Convert the lambda ratio interval onto the fractional-abundance scale:
  #   FA = 100 * lambda_mut / (lambda_ref + lambda_mut)
  #      = 100 * ratio / (1 + ratio)
  # Infinite upper ratio bounds map to 100% fractional abundance.
  tibble(
    fractional_abundance = 100 * mut$lambda / total_lambda,
    ci_low = 100 * ratio_ci$ratio_lower / (1 + ratio_ci$ratio_lower),
    ci_high = ifelse(
      is.infinite(ratio_ci$ratio_upper),
      100,
      100 * ratio_ci$ratio_upper / (1 + ratio_ci$ratio_upper)
    ),
    lambda_ref = ref$lambda,
    lambda_mut = mut$lambda,
    ratio = ratio,
    ratio_low = ratio_ci$ratio_lower,
    ratio_high = ratio_ci$ratio_upper,
    ratio_unbounded = ratio_ci$ratio_unbounded
  )
}

# Select the WT and mutant targets for an assay, requiring different channels.
selected_target_indices <- function(targets, assay) {
  # Pull comparable names and channels out of the raw target metadata.
  names <- vapply(targets, target_name, character(1))
  cleaned <- clean_ddpcr_target(names)
  channels <- vapply(targets, target_channel, integer(1))

  # Identify candidate target entries for the mutant assay and WT reference.
  mut_candidates <- which(cleaned == assay)
  ref_candidates <- which(cleaned == "WT")
  if (length(mut_candidates) == 0L || length(ref_candidates) == 0L) {
    return(NULL)
  }

  # Pair the mutant target with a WT target on the opposite fluorescence channel.
  for (mut_idx in mut_candidates) {
    ref_idx <- ref_candidates[channels[ref_candidates] != channels[mut_idx]][1]
    if (!is.na(ref_idx)) {
      return(list(ref = ref_idx, mut = mut_idx))
    }
  }

  NULL
}

# Extract the first target list from the metadata clusters.
metadata_targets <- function(metadata) {
  # Target metadata is stored inside the first cluster that has at least two targets.
  clusters <- metadata$Clusters %||% list()
  for (cluster in clusters) {
    targets <- cluster$Targets
    if (!is.null(targets) && length(targets) >= 2L) {
      return(targets)
    }
  }
  NULL
}

# Read one well's raw metadata and return analysis count rows.
read_well_count_rows <- function(extract_path, run_row, well_manifest) {
  # Locate this well's peak metadata file within the extracted ddPCR archive.
  assay <- normalise_assay(run_row$assay)
  well <- well_manifest$well[[1]]
  metadata_path <- file.path(extract_path, "PeakMetaData", paste0(well, ".ddmetajson"))

  # The peak metadata is required because it contains per-cluster droplet calls.
  if (!file.exists(metadata_path)) {
    stop("Missing peak metadata: ", metadata_path)
  }

  # Read the JSON metadata and extract the WT/mutant target definitions.
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
  targets <- metadata_targets(metadata)
  if (is.null(targets)) {
    stop("Could not read target metadata for ", run_row$run_id, " ", well)
  }

  # Match this assay's mutant target to its WT reference target.
  selected <- selected_target_indices(targets, assay)
  if (is.null(selected)) {
    stop("Could not select WT and mutant targets for ", run_row$run_id, " ", well)
  }

  # Map the selected WT/mutant targets onto physical channel indices.
  target_names <- vapply(targets, target_name, character(1))
  target_clean <- clean_ddpcr_target(target_names)
  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)
  ch1_idx <- selected_indices[channels[selected_indices] == 1L][1]
  ch2_idx <- selected_indices[channels[selected_indices] == 2L][1]
  if (is.na(ch1_idx) || is.na(ch2_idx)) {
    stop("Could not map selected targets to Ch1 and Ch2 for ", run_row$run_id, " ", well)
  }

  # Initialise the four Ch1/Ch2 partition count categories.
  partition_counts <- c(
    `Ch1+Ch2+` = 0L,
    `Ch1+Ch2-` = 0L,
    `Ch1-Ch2+` = 0L,
    `Ch1-Ch2-` = 0L
  )
  gated_or_unassigned <- 0L

  # Walk through raw clusters and add accepted droplets to partition categories.
  for (cluster in metadata$Clusters %||% list()) {
    droplet_count <- cluster_droplet_count(cluster)
    if (droplet_count == 0L) {
      next
    }

    # Exclude clusters that lack selected target calls or are explicitly unassigned.
    results <- as.character(unlist(cluster$Results, use.names = FALSE))
    if (length(results) < max(selected_indices) || is_true(cluster$Unassigned)) {
      gated_or_unassigned <- gated_or_unassigned + droplet_count
      next
    }

    # Exclude clusters whose selected channels are not simple positive/negative calls.
    ch1 <- results[[ch1_idx]]
    ch2 <- results[[ch2_idx]]
    if (!all(c(ch1, ch2) %in% c("Negative", "Positive"))) {
      gated_or_unassigned <- gated_or_unassigned + droplet_count
      next
    }

    # Convert the Ch1/Ch2 calls into the matching partition-count key.
    key <- paste0(
      ifelse(ch1 == "Positive", "Ch1+", "Ch1-"),
      ifelse(ch2 == "Positive", "Ch2+", "Ch2-")
    )
    partition_counts[[key]] <- partition_counts[[key]] + droplet_count
  }

  # Accepted droplets are those assigned to one of the four partition categories.
  accepted <- sum(partition_counts)
  if (accepted <= 0L) {
    stop("No accepted droplets for active well ", run_row$run_id, " ", well)
  }

  # Convert channel-level partitions into target-index positive counts.
  ch1_positive <- partition_counts[["Ch1+Ch2+"]] + partition_counts[["Ch1+Ch2-"]]
  ch2_positive <- partition_counts[["Ch1+Ch2+"]] + partition_counts[["Ch1-Ch2+"]]
  positives_by_index <- integer(length(targets))
  positives_by_index[ch1_idx] <- ch1_positive
  positives_by_index[ch2_idx] <- ch2_positive

  # Separate the selected target counts into reference and mutant channels.
  ref_positive <- positives_by_index[[selected$ref]]
  mut_positive <- positives_by_index[[selected$mut]]
  ref_negative <- accepted - ref_positive
  mut_negative <- accepted - mut_positive

  # Calculate the official Poisson-corrected fractional abundance metrics.
  fa <- fractional_abundance(
    ref_positive = ref_positive,
    ref_negative = ref_negative,
    mut_positive = mut_positive,
    mut_negative = mut_negative,
    total = accepted
  )

  # Build one row per selected target, preserving raw and derived count fields.
  raw_rows <- tibble(
    target_clean = target_clean[selected_indices],
    target_index = selected_indices,
    raw_target_name = target_names[selected_indices],
    raw_channel = paste0("Ch", channels[selected_indices]),
    raw_role = if_else(selected_indices == selected$mut, "mutant", "reference"),
    `Accepted Droplets` = accepted,
    Positives = positives_by_index[selected_indices],
    Negatives = accepted - positives_by_index[selected_indices],
    `Ch1+Ch2+` = partition_counts[["Ch1+Ch2+"]],
    `Ch1+Ch2-` = partition_counts[["Ch1+Ch2-"]],
    `Ch1-Ch2+` = partition_counts[["Ch1-Ch2+"]],
    `Ch1-Ch2-` = partition_counts[["Ch1-Ch2-"]],
    raw_gated_or_unassigned_droplets = gated_or_unassigned,
    raw_metadata_path = metadata_path
  )

  # Attach manifest fields and keep fractional-abundance metrics on mutant rows.
  raw_rows %>%
    left_join(
      well_manifest %>%
        select(
          run_id, Date = run_date, Well = well, Sample = sample,
          Target = target, ExperimentType = assay, target_clean,
          source_csv, source_ddpcr, source_layout
        ),
      by = "target_clean"
    ) %>%
    mutate(
      `Fractional Abundance` = if_else(
        raw_role == "mutant",
        fa$fractional_abundance,
        NA_real_
      ),
      PoissonFractionalAbundanceMin = if_else(
        raw_role == "mutant",
        fa$ci_low,
        NA_real_
      ),
      PoissonFractionalAbundanceMax = if_else(
        raw_role == "mutant",
        fa$ci_high,
        NA_real_
      ),
      lambda_ref = if_else(
        raw_role == "mutant",
        fa$lambda_ref,
        NA_real_
      ),
      lambda_mut = if_else(
        raw_role == "mutant",
        fa$lambda_mut,
        NA_real_
      ),
      concentration_ratio = if_else(
        raw_role == "mutant",
        fa$ratio,
        NA_real_
      ),
      concentration_ratio_low = if_else(
        raw_role == "mutant",
        fa$ratio_low,
        NA_real_
      ),
      concentration_ratio_high = if_else(
        raw_role == "mutant",
        fa$ratio_high,
        NA_real_
      ),
      concentration_ratio_unbounded = if_else(
        raw_role == "mutant",
        fa$ratio_unbounded,
        NA
      )
    ) %>%
    # Return columns in the order expected by downstream analysis scripts.
    select(
      Sample, Date, Well, Target, ExperimentType,
      `Accepted Droplets`, Positives, Negatives,
      `Ch1+Ch2+`, `Ch1+Ch2-`, `Ch1-Ch2+`, `Ch1-Ch2-`,
      `Fractional Abundance`,
      PoissonFractionalAbundanceMax,
      PoissonFractionalAbundanceMin,
      run_id, target_clean, raw_target_name, raw_channel, raw_role,
      lambda_ref, lambda_mut,
      concentration_ratio, concentration_ratio_low, concentration_ratio_high,
      concentration_ratio_unbounded,
      raw_gated_or_unassigned_droplets, raw_metadata_path,
      source_csv, source_ddpcr, source_layout
  )
}

# Read all active wells for one archived run directory.
read_archive_rows_from_database <- function(raw_root, run_row, sample_manifest) {
  # Resolve the extracted archive directory recorded in the run manifest.
  extract_path <- file.path(raw_root, run_row$archive_contents_relative_dir)
  if (!dir.exists(extract_path)) {
    stop("Missing archive contents directory: ", extract_path)
  }

  # Keep only sample-manifest rows that belong to this run.
  run_manifest <- sample_manifest %>%
    filter(run_id == run_row$run_id) %>%
    mutate(Date = as.Date(run_date))

  # Rebuild all wells from the archive metadata and bind them into one table.
  purrr::map_dfr(
    unique(run_manifest$well),
    # Build and process the manifest rows for this well's WT and mutant targets.
    function(well) {
      well_manifest <- run_manifest %>%
        filter(.data$well == !!well) %>%
        arrange(target_order)

      read_well_count_rows(extract_path, run_row, well_manifest)
    }
  )
}

# Compare raw-derived rows with the active CSV manifest and write validation reports.
validate_raw_rows_against_manifest <- function(raw_rows, sample_manifest, validation_dir) {
  # Ensure the validation report directory exists before writing outputs.
  dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

  # Build the reference table from the active CSV-backed sample manifest.
  reference <- sample_manifest %>%
    transmute(
      run_id,
      Date = as.Date(run_date),
      Well = well,
      Target = target,
      ExperimentType = assay,
      Sample = sample,
      accepted_droplets = as.numeric(accepted_droplets),
      positives = as.numeric(positives),
      negatives = as.numeric(negatives),
      ch1_ch2_pos = as.numeric(ch1_ch2_pos),
      ch1_pos_ch2_neg = as.numeric(ch1_pos_ch2_neg),
      ch1_neg_ch2_pos = as.numeric(ch1_neg_ch2_pos),
      ch1_ch2_neg = as.numeric(ch1_ch2_neg),
      csv_fractional_abundance = as.numeric(fractional_abundance),
      csv_ci_low = as.numeric(poisson_fractional_abundance_min),
      csv_ci_high = as.numeric(poisson_fractional_abundance_max)
    )

  # Build the comparison table from raw-derived archive metadata rows.
  raw_compare <- raw_rows %>%
    transmute(
      run_id,
      Date = as.Date(Date),
      Well,
      Target,
      ExperimentType,
      Sample,
      raw_accepted_droplets = as.numeric(`Accepted Droplets`),
      raw_positives = as.numeric(Positives),
      raw_negatives = as.numeric(Negatives),
      raw_ch1_ch2_pos = as.numeric(`Ch1+Ch2+`),
      raw_ch1_pos_ch2_neg = as.numeric(`Ch1+Ch2-`),
      raw_ch1_neg_ch2_pos = as.numeric(`Ch1-Ch2+`),
      raw_ch1_ch2_neg = as.numeric(`Ch1-Ch2-`),
      calculated_fractional_abundance = `Fractional Abundance`,
      calculated_ci_low = PoissonFractionalAbundanceMin,
      calculated_ci_high = PoissonFractionalAbundanceMax,
      lambda_ref,
      lambda_mut,
      concentration_ratio,
      concentration_ratio_low,
      concentration_ratio_high,
      concentration_ratio_unbounded,
      raw_gated_or_unassigned_droplets
    )

  # Join on row identity so missing and mismatched rows can be detected.
  key_cols <- c("run_id", "Date", "Well", "Target", "ExperimentType", "Sample")
  joined <- full_join(reference, raw_compare, by = key_cols)

  # Calculate count differences and fractional-abundance differences.
  comparison <- joined %>%
    mutate(
      missing_in_manifest = is.na(accepted_droplets),
      missing_in_raw = is.na(raw_accepted_droplets),
      accepted_droplets_diff = raw_accepted_droplets - accepted_droplets,
      positives_diff = raw_positives - positives,
      negatives_diff = raw_negatives - negatives,
      ch1_ch2_pos_diff = raw_ch1_ch2_pos - ch1_ch2_pos,
      ch1_pos_ch2_neg_diff = raw_ch1_pos_ch2_neg - ch1_pos_ch2_neg,
      ch1_neg_ch2_pos_diff = raw_ch1_neg_ch2_pos - ch1_neg_ch2_pos,
      ch1_ch2_neg_diff = raw_ch1_ch2_neg - ch1_ch2_neg,
      csv_vs_calculated_fa_diff = calculated_fractional_abundance - csv_fractional_abundance,
      csv_vs_calculated_ci_low_diff = calculated_ci_low - csv_ci_low,
      csv_vs_calculated_ci_high_diff = calculated_ci_high - csv_ci_high,
      count_difference = missing_in_manifest | missing_in_raw |
        accepted_droplets_diff != 0 |
        positives_diff != 0 |
        negatives_diff != 0 |
        ch1_ch2_pos_diff != 0 |
        ch1_pos_ch2_neg_diff != 0 |
        ch1_neg_ch2_pos_diff != 0 |
        ch1_ch2_neg_diff != 0,
      csv_fa_difference = if_else(
        is.na(csv_vs_calculated_fa_diff),
        FALSE,
        abs(csv_vs_calculated_fa_diff) > 1e-6 |
          abs(csv_vs_calculated_ci_low_diff) > 1e-6 |
          abs(csv_vs_calculated_ci_high_diff) > 1e-6
      )
    )

  # Keep only rows with count differences or CSV-vs-calculated FA differences.
  row_differences <- comparison %>%
    filter(count_difference | csv_fa_difference)

  # Summarise the validation outcome for quick audit.
  summary <- tibble(
    metric = c(
      "manifest_rows",
      "raw_rows",
      "joined_rows",
      "missing_in_manifest",
      "missing_in_raw",
      "rows_with_count_differences",
      "mutant_rows_with_csv_vs_calculated_fractional_abundance_differences",
      "active_wells",
      "raw_gated_or_unassigned_droplets"
    ),
    value = c(
      nrow(reference),
      nrow(raw_compare),
      nrow(joined),
      sum(joined$accepted_droplets %>% is.na()),
      sum(joined$raw_accepted_droplets %>% is.na()),
      sum(row_differences$count_difference, na.rm = TRUE),
      sum(row_differences$csv_fa_difference, na.rm = TRUE),
      n_distinct(raw_rows$run_id, raw_rows$Well),
      sum(raw_rows$raw_gated_or_unassigned_droplets, na.rm = TRUE) / 2
    )
  )

  # Preserve detailed CSV-vs-calculated fractional-abundance comparisons for mutant rows.
  fractional_abundance_comparison <- comparison %>%
    filter(!is.na(calculated_fractional_abundance) | !is.na(csv_fractional_abundance)) %>%
    select(
      all_of(key_cols),
      csv_fractional_abundance,
      csv_ci_low,
      csv_ci_high,
      calculated_fractional_abundance,
      calculated_ci_low,
      calculated_ci_high,
      lambda_ref,
      lambda_mut,
      concentration_ratio,
      concentration_ratio_low,
      concentration_ratio_high,
      concentration_ratio_unbounded,
      csv_vs_calculated_fa_diff,
      csv_vs_calculated_ci_low_diff,
      csv_vs_calculated_ci_high_diff
    )

  # Write validation artefacts beside the analysis outputs.
  readr::write_csv(summary, file.path(validation_dir, "raw_vs_csv_summary.csv"))
  readr::write_csv(row_differences, file.path(validation_dir, "raw_vs_csv_row_differences.csv"))
  readr::write_csv(fractional_abundance_comparison, file.path(validation_dir, "fractional_abundance_method_comparison.csv"))

  # Treat count mismatches as hard failures; FA method differences are reported.
  if (any(row_differences$count_difference, na.rm = TRUE)) {
    stop("Raw-derived droplet counts differ from the active CSV manifest. See raw_vs_csv_row_differences.csv.")
  }

  # Return the summary invisibly so callers can inspect it without console noise.
  invisible(summary)
}

# Load active SNV runs from the raw ddPCR database and optionally validate them.
read_ddpcr_raw_bigdata_from_database <- function(raw_root, validation_dir = NULL) {
  # Locate the run-level and well/target-level manifests.
  runs_path <- file.path(raw_root, "manifests", "runs.csv")
  sample_manifest_path <- file.path(raw_root, "manifests", "sample_manifest.csv")

  # Read only active SNV runs, normalise assay labels, and impose stable order.
  runs <- readr::read_csv(runs_path, show_col_types = FALSE) %>%
    filter(status == "active", experiment == "SNV") %>%
    mutate(
      run_date = as.Date(run_date),
      assay = normalise_assay(assay)
    ) %>%
    arrange(run_date, assay, run_id)

  # Read SNV manifest rows for the supported assays and clean target labels.
  sample_manifest <- readr::read_csv(sample_manifest_path, show_col_types = FALSE) %>%
    filter(experiment == "SNV", assay %in% mutation_list_raw_import) %>%
    mutate(
      run_date = as.Date(run_date),
      assay = normalise_assay(assay),
      target_clean = clean_ddpcr_target(target_clean)
    ) %>%
    arrange(run_date, assay, well, target_order)

  # Rebuild raw count rows from each active run's extracted archive contents.
  raw_rows <- purrr::map_dfr(
    seq_len(nrow(runs)),
    # Process one active run and return its reconstructed well rows.
    function(i) {
      read_archive_rows_from_database(raw_root, runs[i, ], sample_manifest)
    }
  ) %>%
    # Normalise key fields and sort to match downstream expectations.
    mutate(
      Date = as.Date(Date),
      ExperimentType = normalise_assay(ExperimentType)
    ) %>%
    arrange(Date, ExperimentType, Well, Target)

  # Optionally validate raw-derived rows against the active CSV manifest.
  if (!is.null(validation_dir)) {
    validate_raw_rows_against_manifest(raw_rows, sample_manifest, validation_dir)
    readr::write_csv(raw_rows, file.path(validation_dir, "raw_import_bigdata.csv"))
  }

  # Return the reconstructed bigdata table.
  raw_rows
}
