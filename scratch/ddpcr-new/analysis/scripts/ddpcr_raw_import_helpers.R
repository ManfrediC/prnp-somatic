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
  # Remove QuantaSoft probe/channel suffixes from target labels.
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

# QuantaSoft lookup table for rare positive-count confidence bounds.
quantasoft_count_bound_table <- local({
  readr::read_table(
"positive lower95 upper95 lower68 upper68
0 0 2.99600005 0 1.36899996
1 0.0419999994 4.77600002 0.270999998 2.48900008
2 0.303000003 6.40600014 0.867999971 3.84500003
3 0.712000012 7.95100021 1.55900002 5.13600016
4 1.20700002 9.43099976 2.29699993 6.38899994
5 1.75800002 10.8649998 3.06599998 7.61499977
6 2.3499999 12.2639999 3.85700011 8.81999969
7 2.97399998 13.632 4.66300011 10.0120001
8 3.62199998 14.9790001 5.48400021 11.1890001
9 4.29199982 16.3040009 6.31500006 12.3559999
10 4.97900009 17.6140003 7.15500021 13.5150003
11 5.67999983 18.9109993 8.00199986 14.6669998
12 6.39499998 20.1930008 8.85700035 15.8109999
13 7.12099981 21.4650002 9.71700001 16.9510002
14 7.8579998 22.7259998 10.5819998 18.0849991
15 8.60299969 23.9799995 11.4530001 19.2140007
16 9.35599995 25.2259998 12.3269997 20.3390007
17 10.1169996 26.4640007 13.2049999 21.4610004
18 10.8850002 27.6940002 14.0869999 22.5790005
19 11.6590004 28.9190006 14.9720001 23.6930008
20 12.4390001 30.1389999 15.8599997 24.8050003
21 13.224 31.3530006 16.7509995 25.9139996
22 14.0150003 32.5610008 17.6439991 27.0200005
23 14.809 33.7669983 18.5389996 28.125
24 15.6079998 34.9669991 19.4370003 29.2269993
25 16.4120007 36.1619987 20.3369999 30.3269997
26 17.2189999 37.3549995 21.2390003 31.4249992
27 18.0300007 38.5429993 22.1429996 32.5209999
28 18.8439999 39.7290001 23.0480003 33.6160011
29 19.6609993 40.9119987 23.9549999 34.7089996
30 20.4820004 42.0900002 24.8640003 35.7989998
31 21.3059998 43.2659988 25.7740002 36.8889999
32 22.132 44.4389992 26.6860008 37.9770012
33 22.9610004 45.6100006 27.5990009 39.0639992
34 23.7919998 46.7789993 28.5130005 40.1500015
35 24.6259995 47.9449997 29.4290009 41.2340012
36 25.4619999 49.1090012 30.3460007 42.3170013
37 26.3010006 50.269001 31.2639999 43.3989983
38 27.1420002 51.4280014 32.1819992 44.480999
39 27.9839993 52.5859985 33.1020012 45.5610008
40 28.8290005 53.7410011 34.0239983 46.6380005
41 29.6760006 54.8930016 34.9459991 47.7159996
42 30.5240002 56.0449982 35.8689995 48.7929993
43 31.375 57.1940002 36.7929993 49.8689995
44 32.2270012 58.3419991 37.7169991 50.9449997
45 33.0800018 59.4889984 38.6430016 52.019001
46 33.9350014 60.6339989 39.5690002 53.0929985
47 34.7919998 61.776001 40.4970016 54.1650009
48 35.651001 62.9169998 41.4249992 55.2369995
49 36.5099983 64.0579987 42.3540001 56.3079987
50 37.3720016 65.1959991 43.2830009 57.3790016
51 38.2340012 66.3339996 44.2130013 58.4490013
52 39.0979996 67.4700012 45.144001 59.5180016
53 39.9630013 68.6050034 46.0760002 60.5859985
54 40.8289986 69.7389984 47.0079994 61.6539993
55 41.6969986 70.8700027 47.9410019 62.7210007
56 42.5660019 72.0009995 48.8740005 63.7879982
57 43.4360008 73.1309967 49.8079987 64.8539963
58 44.3069992 74.2600021 50.7420006 65.9199982
59 45.1800003 75.387001 51.6769981 66.9850006
60 46.0530014 76.5139999 52.612999 68.0490036
61 46.9269981 77.6399994 53.5489998 69.112999
62 47.8030014 78.7639999 54.4860001 70.1750031
63 48.6790009 79.8880005 55.4230003 71.237999
64 49.5559998 81.0110016 56.3610001 72.3000031
65 50.4339981 82.1330032 57.2989998 73.3619995
66 51.3139992 83.2519989 58.2369995 74.4250031
67 52.1940002 84.3720016 59.1769981 75.4840012
68 53.0750008 85.4909973 60.1160011 76.5449982
69 53.9570007 86.6090012 61.0559998 77.6050034
70 54.8400002 87.7259979 61.9970016 78.6640015
71 55.7229996 88.8430023 62.9370003 79.723999
72 56.6080017 89.9580002 63.8790016 80.7819977
73 57.493 91.072998 64.8199997 81.8410034
74 58.3790016 92.1869965 65.7630005 82.8980026
75 59.2659988 93.3000031 66.7050018 83.9560013
76 60.1529999 94.413002 67.6480026 85.0130005
77 61.0410004 95.5250015 68.5910034 86.0699997
78 61.9300003 96.6360016 69.5350037 87.1259995
79 62.8190002 97.7470016 70.4789963 88.1819992
80 63.7099991 98.8560028 71.4229965 89.237999
81 64.6009979 99.9639969 72.3679962 90.2929993
82 65.4929962 101.071999 73.3130035 91.3479996
83 66.3850021 102.18 74.2580032 92.4029999
84 67.2779999 103.287003 75.2040024 93.4570007
85 68.1719971 104.392998 76.1500015 94.5110016
86 69.0660019 105.499001 77.0960007 95.5650024
87 69.9609985 106.603996 78.0419998 96.6190033
88 70.8560028 107.709 78.9889984 97.6719971
89 71.7519989 108.813004 79.9369965 98.723999
90 72.6490021 109.916 80.8840027 99.7770004
91 73.5459976 111.018997 81.8320007 100.829002
92 74.4440002 112.121002 82.7799988 101.880997
93 75.3420029 113.223 83.7279968 102.932999
94 76.2409973 114.323997 84.677002 103.984001
95 77.1399994 115.425003 85.6259995 105.035004
96 78.0400009 116.525002 86.5749969 106.085999
97 78.9410019 117.624001 87.5250015 107.136002
98 79.8420029 118.723 88.473999 108.186996
99 80.7429962 119.821999 89.4240036 109.237
", show_col_types = FALSE)
})

quantasoft_lambda_cap <- 1000

# Convert positive and negative partition counts to a QuantaSoft-style lambda.
quantasoft_lambda <- function(positive, negative, total) {
  case_when(
    is.na(positive) | is.na(negative) | is.na(total) | total <= 0 ~ NA_real_,
    positive < 0 | negative < 0 | positive > total | negative > total ~ NA_real_,
    negative == 0 ~ quantasoft_lambda_cap,
    negative == total ~ 0,
    TRUE ~ -log(negative / total)
  )
}

# Convert a positive partition count to lambda using the accepted droplet total.
quantasoft_count_to_lambda <- function(positive_count, total) {
  case_when(
    is.na(positive_count) | is.na(total) | total <= 0 ~ NA_real_,
    positive_count <= 0 ~ 0,
    positive_count >= total ~ quantasoft_lambda_cap,
    TRUE ~ -log1p(-positive_count / total)
  )
}

# Reconstruct QuantaSoft positive-count confidence bounds for a partition count.
quantasoft_positive_count_bounds <- function(positive, negative, total, conf.level = 0.95) {
  # Work numerically and initialise bounds for every input row.
  positive <- as.numeric(positive)
  negative <- as.numeric(negative)
  total <- as.numeric(total)
  lower <- upper <- rep(NA_real_, length(positive))
  ok <- !is.na(positive) & !is.na(negative) & !is.na(total) & total > 0

  # Select the QuantaSoft lookup columns and normal approximation multiplier.
  if (abs(conf.level - 0.95) < 1e-8) {
    lower_column <- "lower95"
    upper_column <- "upper95"
    z <- 1.96
  } else if (abs(conf.level - 0.68) < 1e-8) {
    lower_column <- "lower68"
    upper_column <- "upper68"
    z <- 1
  } else {
    stop("QuantaSoft reconstruction supports only 0.95 and 0.68 confidence levels.")
  }

  # For fewer than 100 positives, use the rare-count lookup table directly.
  rare_positive <- ok & positive < 100
  if (any(rare_positive)) {
    idx <- as.integer(positive[rare_positive]) + 1L
    lower[rare_positive] <- quantasoft_count_bound_table[[lower_column]][idx]
    upper[rare_positive] <- quantasoft_count_bound_table[[upper_column]][idx]
  }

  # For near-saturated wells, mirror the rare-count lookup around negatives.
  rare_negative <- ok & !rare_positive & negative < 100
  if (any(rare_negative)) {
    idx <- as.integer(negative[rare_negative]) + 1L
    lower[rare_negative] <- total[rare_negative] - quantasoft_count_bound_table[[upper_column]][idx]
    upper[rare_negative] <- total[rare_negative] - quantasoft_count_bound_table[[lower_column]][idx]
  }

  # For the middle region, use QuantaSoft's normal approximation.
  normal_region <- ok & !rare_positive & !rare_negative
  if (any(normal_region)) {
    sigma <- sqrt(positive[normal_region] * negative[normal_region] / total[normal_region])
    lower[normal_region] <- positive[normal_region] - z * sigma
    upper[normal_region] <- positive[normal_region] + z * sigma
  }

  # Clamp bounds back to the physically possible count range.
  tibble(
    lower = pmax(0, pmin(total, lower)),
    upper = pmax(0, pmin(total, upper))
  )
}

# Convert QuantaSoft-style count bounds into lambda bounds.
quanta_like_lambda_bounds <- function(positive, negative, total, conf.level = 0.95) {
  # First reconstruct QuantaSoft-like bounds on the count scale.
  count_bounds <- quantasoft_positive_count_bounds(positive, negative, total, conf.level = conf.level)

  # Then transform those count bounds to the Poisson lambda scale.
  tibble(
    lower = quantasoft_count_to_lambda(count_bounds$lower, total),
    upper = quantasoft_count_to_lambda(count_bounds$upper, total)
  )
}

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
    return(quantasoft_lambda_cap)
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
      quantasoft_lambda_cap,
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

# Recreate QuantaSoft-like fractional abundance and propagated confidence bounds.
quanta_like_fractional_abundance <- function(ref_positive, ref_negative, mut_positive, mut_negative, total) {
  # Recalculate QuantaSoft-style lambdas and the fractional abundance point estimate.
  ref_lambda <- quantasoft_lambda(ref_positive, ref_negative, total)
  mut_lambda <- quantasoft_lambda(mut_positive, mut_negative, total)
  fa <- fractional_abundance_from_lambdas(ref_lambda, mut_lambda)

  # Preserve output shape when the point estimate is unavailable.
  if (is.na(fa) || total <= 0) {
    return(tibble(
      quanta_like_fractional_abundance = NA_real_,
      quanta_like_ci_low = NA_real_,
      quanta_like_ci_high = NA_real_
    ))
  }

  # Reconstruct QuantaSoft-like lambda bounds for each channel.
  ref_bounds <- quanta_like_lambda_bounds(ref_positive, ref_negative, total)
  mut_bounds <- quanta_like_lambda_bounds(mut_positive, mut_negative, total)

  # Handle the all-zero concentration case explicitly.
  denominator <- ref_lambda + mut_lambda
  if (denominator <= 0) {
    return(tibble(
      quanta_like_fractional_abundance = 0,
      quanta_like_ci_low = 0,
      quanta_like_ci_high = 0
    ))
  }

  # Propagate reference and mutant lambda uncertainty to fractional abundance.
  ref_half_width <- (ref_bounds$upper - ref_bounds$lower) / 2
  mut_half_width <- (mut_bounds$upper - mut_bounds$lower) / 2
  d_ref <- -mut_lambda / denominator^2
  d_mut <- ref_lambda / denominator^2
  error_fraction <- sqrt((d_ref * ref_half_width)^2 + (d_mut * mut_half_width)^2)

  # Return a lower bound clipped at zero and an unclipped upper bound.
  tibble(
    quanta_like_fractional_abundance = fa,
    quanta_like_ci_low = pmax(0, fa - 100 * error_fraction),
    quanta_like_ci_high = fa + 100 * error_fraction
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

# Read one well's raw metadata and return QuantaSoft-shaped count rows.
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

  # Calculate both analysis and QuantaSoft-like fractional abundance metrics.
  fa <- fractional_abundance(
    ref_positive = ref_positive,
    ref_negative = ref_negative,
    mut_positive = mut_positive,
    mut_negative = mut_negative,
    total = accepted
  )
  quanta_fa <- quanta_like_fractional_abundance(
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
      ),
      quanta_like_fractional_abundance = if_else(
        raw_role == "mutant",
        quanta_fa$quanta_like_fractional_abundance,
        NA_real_
      ),
      quanta_like_ci_low = if_else(
        raw_role == "mutant",
        quanta_fa$quanta_like_ci_low,
        NA_real_
      ),
      quanta_like_ci_high = if_else(
        raw_role == "mutant",
        quanta_fa$quanta_like_ci_high,
        NA_real_
      )
    ) %>%
    # Return columns in the QuantaSoft-shaped order expected downstream.
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
      quanta_like_fractional_abundance, quanta_like_ci_low, quanta_like_ci_high,
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
      quanta_like_fractional_abundance,
      quanta_like_ci_low,
      quanta_like_ci_high,
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

  # Preserve detailed fractional-abundance method comparisons for mutant rows.
  fa_method_comparison <- comparison %>%
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
      quanta_like_fractional_abundance,
      quanta_like_ci_low,
      quanta_like_ci_high,
      csv_vs_calculated_fa_diff,
      csv_vs_calculated_ci_low_diff,
      csv_vs_calculated_ci_high_diff
    )

  # Write validation artefacts beside the analysis outputs.
  readr::write_csv(summary, file.path(validation_dir, "raw_vs_csv_summary.csv"))
  readr::write_csv(row_differences, file.path(validation_dir, "raw_vs_csv_row_differences.csv"))
  readr::write_csv(fa_method_comparison, file.path(validation_dir, "fa_method_comparison.csv"))

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

  # Return the reconstructed QuantaSoft-shaped bigdata table.
  raw_rows
}
