#!/usr/bin/env Rscript

# Build compact DNA-quality sample evidence table.
# Run from the repository root.

# -------------------------------------------------------------------------
# Fixed input/output paths
# -------------------------------------------------------------------------

MANIFEST_PATH <- "raw/dna_quality/metadata/dna_quality_manifest.tsv"
DDPCR_XLSX <- "results/ddPCR/SNV_data_final.xlsx"
SEQUENCING_TSV <- "results/sequencing_qc/sequencing_metrics_per_sample.tsv"
OUTPUT_PATH <- "results/dna_quality/sample_quality_evidence_table.tsv"
LATEX_OUTPUT_PATH <- "manuscript/tables/dna_quality/sample_quality_evidence_table.tex"

# -------------------------------------------------------------------------
# Basic readers
# -------------------------------------------------------------------------

read_tsv <- function(path) {
  # Read UTF-8/BOM-aware TSV files while preserving original column names.
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fileEncoding = "UTF-8-BOM"
  )
}

read_csv_export <- function(path) {
  # Read TapeStation CSV exports while preserving vendor column names and units.
  # The exports use a Windows/Latin-1 micro sign in unit labels such as ng/ul,
  # so read the lines with latin1 encoding before parsing the CSV text.
  lines <- readLines(file(path, encoding = "latin1"), warn = FALSE)

  # Parse from text to avoid noisy "incomplete final line" warnings in exports.
  read.csv(
    text = paste(lines, collapse = "\n"),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

read_xlsx <- function(path) {
  # Keep Excel parsing explicit so missing dependencies fail with a useful message.
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required. Install it with install.packages('readxl').")
  }

  # Return a plain data.frame so downstream code does not depend on tibble behavior.
  as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------
# Small helpers
# -------------------------------------------------------------------------

find_column <- function(df, patterns, required = FALSE) {
  # Search case-insensitively while returning the original exported column name.
  column_names <- tolower(names(df))

  # Try each acceptable pattern in priority order and return the first match.
  for (pattern in patterns) {
    hit <- grep(pattern, column_names, perl = TRUE)
    if (length(hit) > 0) {
      return(names(df)[hit[1]])
    }
  }

  # Stop only for columns that are impossible to continue without.
  if (required) {
    stop("Missing expected column: ", paste(patterns, collapse = " | "))
  }

  NA_character_
}

numeric_part <- function(x) {
  # Normalize common empty/missing values before numeric conversion.
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN")] <- NA_character_

  # Remove formatting characters and inequality qualifiers, keeping the number.
  x <- gsub(",", "", x, fixed = TRUE)
  x <- sub("^[<>]", "", x)

  # Suppress warnings because missing/non-numeric exports should become NA.
  suppressWarnings(as.numeric(x))
}

reported_number <- function(x, scale = 1) {
  # Preserve the original inequality qualifier when values are out of range.
  raw <- trimws(as.character(x))
  raw[raw %in% c("", "NA", "NaN")] <- NA_character_

  # Capture ">" or "<" so values such as >60000 remain interpretable.
  qualifier <- ifelse(
    !is.na(raw) & grepl("^[<>]", raw),
    substr(raw, 1, 1),
    ""
  )

  # Convert to a scaled numeric value, for example pg/ul to ng/ul.
  value <- numeric_part(raw) * scale

  # Emit compact text for the manuscript-facing evidence table.
  ifelse(
    is.na(value),
    NA_character_,
    paste0(qualifier, format(signif(value, 6), scientific = FALSE, trim = TRUE))
  )
}

normalize_well <- function(x) {
  # Normalize well labels such as A01 and A1 to the same value for matching.
  x <- toupper(trimws(as.character(x)))

  # Keep non-standard labels unchanged so unexpected exports remain visible.
  ifelse(
    grepl("^[A-Z]+[0-9]+$", x),
    paste0(
      sub("^([A-Z]+).*", "\\1", x),
      as.integer(sub("^[A-Z]+0*", "", x))
    ),
    x
  )
}

control_name_for_ddpcr <- function(sample_id) {
  # ddPCR final results use Control1-style names, while the manifest uses Ctrl1.
  sub("^Ctrl", "Control", sample_id)
}

# -------------------------------------------------------------------------
# Manifest selection
# -------------------------------------------------------------------------

select_manifest_row <- function(manifest, sample_id, stage, require_analysis_use = FALSE) {
  # Select candidate manifest rows for exactly one canonical sample and stage.
  rows <- manifest[
    manifest$sample_id == sample_id & manifest$stage == stage,
    ,
    drop = FALSE
  ]

  # Missing rows stay explicit as NULL so callers can write blank output fields.
  if (nrow(rows) == 0) {
    return(NULL)
  }

  # For pre-capture evidence, prefer rows approved for default analysis use.
  if (require_analysis_use && "analysis_use" %in% names(rows)) {
    preferred <- rows[tolower(rows$analysis_use) == "yes", , drop = FALSE]

    if (nrow(preferred) > 0) {
      rows <- preferred
    }
  }

  # The reviewed manifest should be unique here; choose the first row defensively.
  rows[1, , drop = FALSE]
}

match_tapestation_row <- function(df, manifest_row) {
  # Locate the two most stable identifiers available in TapeStation exports.
  sample_col <- find_column(df, c("^sample description$", "sample.*description"))
  well_col <- find_column(df, c("^well$", "^wellid$", "^well id$"))

  # Start with all rows, then narrow candidates using sample description and well.
  candidates <- seq_len(nrow(df))

  # Prefer exact sample-description matching because it is closest to the manifest.
  if (!is.na(sample_col)) {
    sample_label <- trimws(manifest_row$tapestation_sample_description)
    exact <- candidates[trimws(as.character(df[[sample_col]])) == sample_label]

    # Some exports use canonical sample IDs instead of the raw Tapestation label.
    if (length(exact) == 0) {
      exact <- candidates[trimws(as.character(df[[sample_col]])) == manifest_row$sample_id]
    }

    # Keep all exact matches for now; the well filter below can disambiguate.
    if (length(exact) > 0) {
      candidates <- exact
    }
  }

  # If sample labels are not unique, use normalized well labels as a second key.
  if (!is.na(well_col) && length(candidates) > 1) {
    well_label <- normalize_well(manifest_row$tapestation_well)
    well_match <- candidates[normalize_well(df[[well_col]][candidates]) == well_label]

    if (length(well_match) > 0) {
      candidates <- well_match
    }
  }

  # No candidate means this manifest row cannot be linked to this export table.
  if (length(candidates) == 0) {
    return(NULL)
  }

  # Return exactly one matched export row for metric extraction.
  df[candidates[1], , drop = FALSE]
}

# -------------------------------------------------------------------------
# TapeStation extraction
# -------------------------------------------------------------------------

extract_sample_concentration_ng_ul <- function(path, manifest_row) {
  # Missing source files become blank evidence fields rather than hard failures.
  if (is.na(path) || !file.exists(path)) {
    return(NA_character_)
  }

  # Read the sample table and link it back to the reviewed manifest row.
  df <- read_csv_export(path)
  row <- match_tapestation_row(df, manifest_row)

  # If the row cannot be matched, leave the concentration blank for review.
  if (is.null(row)) {
    return(NA_character_)
  }

  # Locate the vendor concentration column regardless of D1000/HSD1000 unit text.
  conc_col <- find_column(row, c("^conc\\.", "conc.*\\[.*ng", "conc.*\\[.*pg", "conc"))

  # No concentration column means this table cannot support the metric.
  if (is.na(conc_col)) {
    return(NA_character_)
  }

  # Normalize HSD1000 pg/ul values to ng/ul for a single final-table unit.
  scale <- if (grepl("pg", tolower(conc_col), fixed = TRUE)) {
    1 / 1000
  } else {
    1
  }

  # Return a compact reported value, preserving any out-of-range qualifier.
  reported_number(row[[conc_col]], scale = scale)
}

extract_dominant_peak_bp <- function(path, manifest_row) {
  # Missing peak tables leave dominant peak blank rather than dropping the sample.
  if (is.na(path) || !file.exists(path)) {
    return(NA_character_)
  }

  # Peak tables can have multiple rows per sample, including markers and peaks.
  df <- read_csv_export(path)

  # Find sample and comment columns so we can isolate true sample peaks.
  sample_col <- find_column(df, c("^sample description$", "sample.*description"))
  comment_col <- find_column(df, c("peak comment", "observations"))

  # Restrict to the reviewed sample label from the manifest.
  if (!is.na(sample_col)) {
    df <- df[
      trimws(as.character(df[[sample_col]])) ==
        trimws(manifest_row$tapestation_sample_description),
      ,
      drop = FALSE
    ]
  }

  # If no peak row matches the sample, the evidence field stays blank.
  if (nrow(df) == 0) {
    return(NA_character_)
  }

  # Remove lower/upper marker rows so marker peaks are not treated as DNA peaks.
  if (!is.na(comment_col)) {
    comments <- tolower(as.character(df[[comment_col]]))
    df <- df[!grepl("marker", comments, fixed = TRUE), , drop = FALSE]
  }

  # If only marker rows existed, there is no dominant sample peak to report.
  if (nrow(df) == 0) {
    return(NA_character_)
  }

  # Use exported bp size as the calibrated fragment-size value.
  size_col <- find_column(df, c("size.*bp"), required = TRUE)

  # Rank peaks by integrated area when available, otherwise by concentration.
  area_col <- find_column(df, c("%.*integrated.*area", "integrated.*area"))
  conc_col <- find_column(df, c("calibrated conc", "assigned conc", "conc"))

  # Compute the ranking value from the best available peak-strength field.
  if (!is.na(area_col)) {
    rank_value <- numeric_part(df[[area_col]])
  } else if (!is.na(conc_col)) {
    rank_value <- numeric_part(df[[conc_col]])
  } else {
    rank_value <- rep(NA_real_, nrow(df))
  }

  # If no ranking metric exists, fall back to the first non-marker peak row.
  chosen <- if (all(is.na(rank_value))) {
    1
  } else {
    which.max(replace(rank_value, is.na(rank_value), -Inf))
  }

  # Return the dominant peak size as reported by the TapeStation export.
  reported_number(df[[size_col]][chosen])
}

extract_average_region_bp <- function(path, manifest_row) {
  # Missing region tables leave average region size blank rather than failing.
  if (is.na(path) || !file.exists(path)) {
    return(NA_character_)
  }

  # Read the compact region table and link it back to the reviewed manifest row.
  df <- read_csv_export(path)
  row <- match_tapestation_row(df, manifest_row)

  # If no region row matches the sample, the metric remains blank.
  if (is.null(row)) {
    return(NA_character_)
  }

  # Use the vendor-reported average bp size across the defined region.
  avg_col <- find_column(row, c("average size.*bp", "average.*bp"))

  # Some tables may not export average region size.
  if (is.na(avg_col)) {
    return(NA_character_)
  }

  # Preserve any exported out-of-range qualifier in the final table.
  reported_number(row[[avg_col]])
}

build_tapestation_summary <- function(manifest, sample_ids) {
  # Build one TapeStation evidence row per canonical manifest sample.
  rows <- lapply(sample_ids, function(sample_id) {
    # Pre-capture rows are the default extraction/decontamination evidence.
    pre <- select_manifest_row(
      manifest,
      sample_id,
      stage = "pre_capture",
      require_analysis_use = TRUE
    )

    # Submission QC rows provide downstream sample-prep context.
    submission <- select_manifest_row(
      manifest,
      sample_id,
      stage = "submission_qc",
      require_analysis_use = FALSE
    )

    # Extract the compact TapeStation block selected for the final evidence table.
    data.frame(
      sample_id = sample_id,
      batch_id = if (is.null(pre)) NA_character_ else pre$batch_id,
      pre_capture_tapestation_concentration_ng_ul =
        if (is.null(pre)) NA_character_
        else extract_sample_concentration_ng_ul(pre$sample_table_path, pre),
      pre_capture_dominant_peak_bp =
        if (is.null(pre)) NA_character_
        else extract_dominant_peak_bp(pre$peak_table_path, pre),
      pre_capture_average_region_bp =
        if (is.null(pre)) NA_character_
        else extract_average_region_bp(pre$region_table_path, pre),
      submission_qc_tapestation_concentration_ng_ul =
        if (is.null(submission)) NA_character_
        else extract_sample_concentration_ng_ul(submission$sample_table_path, submission),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# ddPCR summary
# -------------------------------------------------------------------------

build_ddpcr_summary <- function(sample_ids) {
  # If the curated ddPCR workbook is absent, keep all samples visible for review.
  if (!file.exists(DDPCR_XLSX)) {
    return(data.frame(
      sample_id = sample_ids,
      ddpcr_max_fr_accepted_droplets = NA_real_,
      ddpcr_amplification_review_status = "review",
      stringsAsFactors = FALSE
    ))
  }

  # Read the curated ddPCR table that contains n_total_droplets.
  ddpcr <- read_xlsx(DDPCR_XLSX)
  names(ddpcr) <- trimws(names(ddpcr))

  # These columns are required to apply the frontal accepted-droplet rule.
  required <- c("participant", "brain_region", "n_total_droplets")
  missing <- setdiff(required, names(ddpcr))

  # Stop on schema drift so we do not silently report incorrect evidence.
  if (length(missing) > 0) {
    stop("Missing required ddPCR columns: ", paste(missing, collapse = ", "))
  }

  # Summarize frontal cortex accepted droplets for each manifest sample.
  rows <- lapply(sample_ids, function(sample_id) {
    # Harmonize control naming between the manifest and ddPCR workbook.
    participant <- control_name_for_ddpcr(sample_id)

    # The TapeStation/SureSelect samples are frontal cortex, so use fr ddPCR rows.
    sample_rows <- ddpcr[
      ddpcr$participant == participant &
        tolower(trimws(ddpcr$brain_region)) == "fr",
      ,
      drop = FALSE
    ]

    # n_total_droplets is the accepted-droplet count used for amplification fitness.
    droplets <- numeric_part(sample_rows$n_total_droplets)

    # One qualifying frontal assay is enough, so use the maximum accepted droplets.
    max_droplets <- if (length(droplets) == 0 || all(is.na(droplets))) {
      NA_real_
    } else {
      max(droplets, na.rm = TRUE)
    }

    # Mark missing or sub-threshold linkage for review rather than dropping rows.
    data.frame(
      sample_id = sample_id,
      ddpcr_max_fr_accepted_droplets = max_droplets,
      ddpcr_amplification_review_status =
        if (!is.na(max_droplets) && max_droplets >= 10000) "pass" else "review",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# Sequencing summary
# -------------------------------------------------------------------------

build_sequencing_summary <- function(sample_ids) {
  # If sequencing metrics are absent, keep all samples visible for review.
  if (!file.exists(SEQUENCING_TSV)) {
    return(data.frame(
      sample_id = sample_ids,
      seq_coding_pct_bases_ge_100x = NA_real_,
      seq_coding_mean_depth = NA_real_,
      seq_coding_p20_depth = NA_real_,
      sequencing_review_status = "review",
      stringsAsFactors = FALSE
    ))
  }

  # Read the canonical upstream sequencing metrics table.
  sequencing <- read_tsv(SEQUENCING_TSV)

  # These columns define PRNP coding-region coverage at the configured threshold.
  required <- c(
    "sample_id",
    "cov_thr_x1",
    "coding_pct_bases_ge_x1",
    "coding_mean_depth",
    "coding_p20_depth"
  )

  # Stop on schema drift so coverage values are not misinterpreted.
  missing <- setdiff(required, names(sequencing))

  # Missing required columns means the upstream sequencing QC table changed.
  if (length(missing) > 0) {
    stop("Missing required sequencing columns: ", paste(missing, collapse = ", "))
  }

  # Summarize PRNP coding coverage for each manifest sample.
  rows <- lapply(sample_ids, function(sample_id) {
    # Select rows for this exact canonical sample ID.
    sample_rows <- sequencing[sequencing$sample_id == sample_id, , drop = FALSE]

    # Keep only rows where the first configured coverage threshold is 100x.
    sample_rows <- sample_rows[
      numeric_part(sample_rows$cov_thr_x1) == 100,
      ,
      drop = FALSE
    ]

    # Missing sequencing rows become review status in the final evidence table.
    if (nrow(sample_rows) == 0) {
      return(data.frame(
        sample_id = sample_id,
        seq_coding_pct_bases_ge_100x = NA_real_,
        seq_coding_mean_depth = NA_real_,
        seq_coding_p20_depth = NA_real_,
        sequencing_review_status = "review",
        stringsAsFactors = FALSE
      ))
    }

    # The upstream table should be unique per sample; choose the first row defensively.
    row <- sample_rows[1, , drop = FALSE]

    # Carry coverage as continuous downstream context, not as a pass/fail gate.
    data.frame(
      sample_id = sample_id,
      seq_coding_pct_bases_ge_100x = numeric_part(row$coding_pct_bases_ge_x1),
      seq_coding_mean_depth = numeric_part(row$coding_mean_depth),
      seq_coding_p20_depth = numeric_part(row$coding_p20_depth),
      sequencing_review_status = "available",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# Final status and notes
# -------------------------------------------------------------------------

combine_downstream_assays_worked <- function(ddpcr_status, sequencing_status) {
  # A sample has functional downstream evidence when ddPCR passes and sequencing exists.
  ifelse(
    ddpcr_status == "pass" & sequencing_status == "available",
    "yes",
    "review"
  )
}

compose_notes <- function(row) {
  # Accumulate short review notes without hiding which evidence layer is missing.
  notes <- character()

  # ddPCR review usually means missing frontal linkage or accepted droplets <10000.
  if (row$ddpcr_amplification_review_status == "review") {
    notes <- c(notes, "ddPCR frontal accepted-droplet linkage needs review")
  }

  # Sequencing review usually means missing PRNP coding coverage linkage.
  if (row$sequencing_review_status == "review") {
    notes <- c(notes, "sequencing PRNP coding coverage linkage needs review")
  }

  # Blank notes are intentional for rows with complete downstream linkage.
  paste(notes, collapse = "; ")
}

# -------------------------------------------------------------------------
# Escape values for safe use in a raw LaTeX table.
# -------------------------------------------------------------------------

latex_escape <- function(x) {
  # Convert values to text and keep missing values visually blank.
  x <- as.character(x)
  x[is.na(x)] <- ""

  # Escape the small set of LaTeX-special characters expected in table cells.
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)

  x
}

# -------------------------------------------------------------------------
# Order manuscript samples as CJD1, CJD2, ..., then Control1, Control2, ...
# -------------------------------------------------------------------------

sample_sort_key <- function(sample_id) {
  # Put CJD samples first, sorted by numeric suffix.
  group_rank <- ifelse(grepl("^CJD[0-9]+$", sample_id), 1, 2)

  # Sort controls after CJD samples, also by numeric suffix.
  number <- suppressWarnings(as.integer(sub("^[A-Za-z]+", "", sample_id)))
  number[is.na(number)] <- 999999L

  # Include the original ID as a final stable tie-breaker.
  order(group_rank, number, sample_id)
}

# -------------------------------------------------------------------------
# Convert manifest control IDs to manuscript display names.
# -------------------------------------------------------------------------

sample_display_name <- function(sample_id) {
  # The manuscript table should use Control1-style labels.
  sub("^Ctrl", "Control", sample_id)
}

# -------------------------------------------------------------------------
# Build the concise manuscript-facing table from the full evidence table.
# -------------------------------------------------------------------------

build_latex_table_data <- function(evidence) {
  # Keep only columns that are useful for a compact manuscript table.
  table <- evidence[, c(
    "sample_id",
    "pre_capture_tapestation_concentration_ng_ul",
    "pre_capture_dominant_peak_bp",
    "pre_capture_average_region_bp",
    "submission_qc_tapestation_concentration_ng_ul",
    "ddpcr_max_fr_accepted_droplets",
    "seq_coding_pct_bases_ge_100x",
    "seq_coding_mean_depth",
    "downstream_assays_worked"
  ), drop = FALSE]

  # Apply the requested manuscript sample order before changing display names.
  table <- table[sample_sort_key(table$sample_id), , drop = FALSE]

  # Display controls as Control1, Control2, ... in the manuscript table.
  table$sample_id <- sample_display_name(table$sample_id)

  # Round submission-QC concentration to one decimal place.
  table$submission_qc_tapestation_concentration_ng_ul <- ifelse(
    is.na(table$submission_qc_tapestation_concentration_ng_ul),
    NA,
    sprintf(
      "%.1f",
      round(as.numeric(table$submission_qc_tapestation_concentration_ng_ul), 1)
    )
  )

  # Round PRNP coding bases covered at >=100x to one decimal place.
  table$seq_coding_pct_bases_ge_100x <- ifelse(
    is.na(table$seq_coding_pct_bases_ge_100x),
    NA,
    sprintf("%.1f", round(as.numeric(table$seq_coding_pct_bases_ge_100x), 1))
  )

  # Round mean PRNP coding-region coverage to the nearest integer x-fold value.
  table$seq_coding_mean_depth <- ifelse(
    is.na(table$seq_coding_mean_depth),
    NA,
    as.character(round(as.numeric(table$seq_coding_mean_depth)))
  )

  # Use short, intuitive headers that are suitable for manuscript polishing later.
  names(table) <- c(
    "Sample",
    "Pre-capture conc. (ng/ul)",
    "Main fragment (bp)",
    "Avg. fragment region (bp)",
    "Submission conc. (ng/ul)",
    "ddPCR droplets",
    "PRNP coding bases >=100x (%)",
    "Mean PRNP coverage (x)",
    "ddPCR + sequencing"
  )

  table
}

# -------------------------------------------------------------------------
# Write a minimal raw LaTeX table that can be styled later.
# -------------------------------------------------------------------------

write_latex_table <- function(table, path) {
  # Ensure the manuscript table output directory exists.
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  # Escape headers and cells before writing raw LaTeX.
  escaped <- as.data.frame(lapply(table, latex_escape), stringsAsFactors = FALSE)
  headers <- latex_escape(names(table))

  # Use a simple all-left-aligned tabular for now; manuscript styling can follow.
  align <- paste(rep("l", ncol(escaped)), collapse = "")
  lines <- c(
    "% Auto-generated by src/dna_quality/build_sample_quality_evidence_table.R",
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste(headers, collapse = " & "),
    "\\\\",
    "\\hline"
  )

  # Add one LaTeX row per sample.
  for (i in seq_len(nrow(escaped))) {
    lines <- c(lines, paste(as.character(escaped[i, ]), collapse = " & "), "\\\\")
  }

  # Close the raw table environment.
  lines <- c(lines, "\\hline", "\\end{tabular}", "")

  writeLines(lines, path, useBytes = TRUE)
}

# -------------------------------------------------------------------------
# Main workflow
# -------------------------------------------------------------------------

# Load the reviewed DNA-quality manifest that defines canonical samples and paths.
manifest <- read_tsv(MANIFEST_PATH)

# Preserve one final evidence row per canonical manifest sample.
sample_ids <- sort(unique(manifest$sample_id))

# Carry sample identity and group directly from the reviewed manifest.
identity <- unique(manifest[, c("sample_id", "sample_group"), drop = FALSE])
identity <- identity[match(sample_ids, identity$sample_id), , drop = FALSE]

# Build each evidence layer independently before joining.
tapestation <- build_tapestation_summary(manifest, sample_ids)
ddpcr <- build_ddpcr_summary(sample_ids)
sequencing <- build_sequencing_summary(sample_ids)

# Join evidence layers by canonical sample ID while keeping all manifest samples.
evidence <- Reduce(
  function(x, y) merge(x, y, by = "sample_id", all.x = TRUE, sort = FALSE),
  list(identity, tapestation, ddpcr, sequencing)
)

# Add one compact downstream functional-evidence status for quick review.
evidence$downstream_assays_worked <- combine_downstream_assays_worked(
  evidence$ddpcr_amplification_review_status,
  evidence$sequencing_review_status
)

# Keep the output schema manuscript-friendly.
final_columns <- c(
  "sample_id",
  "sample_group",
  "batch_id",
  "pre_capture_tapestation_concentration_ng_ul",
  "pre_capture_dominant_peak_bp",
  "pre_capture_average_region_bp",
  "submission_qc_tapestation_concentration_ng_ul",
  "ddpcr_max_fr_accepted_droplets",
  "ddpcr_amplification_review_status",
  "seq_coding_pct_bases_ge_100x",
  "seq_coding_mean_depth",
  "seq_coding_p20_depth",
  "sequencing_review_status",
  "downstream_assays_worked"
)

# Ensure the canonical output directory exists.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

# Write a clean TSV that can be reviewed, plotted, or imported into manuscript code.
write.table(
  evidence[, final_columns, drop = FALSE],
  file = OUTPUT_PATH,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# Print the output path for interactive command-line use.
message("Wrote: ", OUTPUT_PATH)

# Build and write the raw manuscript LaTeX table.
latex_table <- build_latex_table_data(evidence)
write_latex_table(latex_table, LATEX_OUTPUT_PATH)
message("Wrote: ", LATEX_OUTPUT_PATH)
