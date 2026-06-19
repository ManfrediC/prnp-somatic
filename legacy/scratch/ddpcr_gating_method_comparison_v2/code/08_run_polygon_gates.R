source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

setup_v2_dirs()

method_id <- "polygon_control_hull_gates"
method_variant <- "control-derived convex-hull polygon gates"
class_names <- c("NN", "WT", "MUT", "DP")
min_hull_points <- 30L
polygon_expansion <- 1.03
ellipse_fallback_points <- 48L
plot_sample_per_group <- 400L
plot_sample_per_well <- 120L
random_seed <- 20260606L

set.seed(random_seed)

json_text <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE))
}

parameter_json <- json_text(list(
  method = "control_convex_hull_polygon_gates",
  control_source = "NTC_WT_control_positive_control_json_classes",
  trim = "control_geometry_mahalanobis_radius",
  min_hull_points = min_hull_points,
  polygon_expansion = polygon_expansion,
  sparse_class_fallback = "control_geometry_ellipse_polygon",
  random_seed = random_seed
))

expand_polygon <- function(vertices, centre_row, expansion) {
  vertices %>%
    mutate(
      ch1_amplitude = centre_row$centre_ch1 + (ch1_amplitude - centre_row$centre_ch1) * expansion,
      ch2_amplitude = centre_row$centre_ch2 + (ch2_amplitude - centre_row$centre_ch2) * expansion
    )
}

ellipse_fallback_polygon <- function(geometry, assay_name, class_name) {
  ellipse_points_from_geometry(
    centres = geometry$centres %>% filter(assay == assay_name, target_class == class_name),
    covariance_rows = geometry$covariance_rows %>% filter(assay == assay_name),
    gate_radii = geometry$gate_radii %>% filter(assay == assay_name),
    points = ellipse_fallback_points
  ) %>%
    select(ch1_amplitude, ch2_amplitude)
}

build_class_polygon <- function(control_droplets, geometry, assay_name, class_name) {
  centre_row <- geometry$centres %>%
    filter(assay == assay_name, target_class == class_name) %>%
    slice(1)
  cov_mat <- covariance_matrix_for(geometry$covariance_rows, assay_name, class_name)
  radius <- geometry$gate_radii %>%
    filter(assay == assay_name, target_class == class_name) %>%
    pull(radius_d2) %>%
    first()

  class_points <- control_droplets %>%
    filter(assay == assay_name, target_class == class_name) %>%
    select(ch1_amplitude, ch2_amplitude) %>%
    filter(if_all(everything(), is.finite)) %>%
    distinct()

  source <- "control_geometry_ellipse_polygon"
  trimmed_n <- 0L
  vertices <- ellipse_fallback_polygon(geometry, assay_name, class_name)

  if (nrow(class_points) >= min_hull_points) {
    d2 <- mahalanobis(
      as.matrix(class_points),
      center = c(centre_row$centre_ch1, centre_row$centre_ch2),
      cov = cov_mat
    )
    trimmed <- class_points[is.finite(d2) & d2 <= radius, , drop = FALSE]
    trimmed_n <- nrow(trimmed)

    if (nrow(trimmed) >= 3L &&
        sd(trimmed$ch1_amplitude, na.rm = TRUE) > 0 &&
        sd(trimmed$ch2_amplitude, na.rm = TRUE) > 0) {
      hull_idx <- chull(trimmed$ch1_amplitude, trimmed$ch2_amplitude)
      if (length(hull_idx) >= 3L) {
        vertices <- trimmed[hull_idx, , drop = FALSE]
        source <- "control_convex_hull_mahalanobis_trimmed"
      }
    }
  }

  vertices %>%
    expand_polygon(centre_row, polygon_expansion) %>%
    mutate(
      assay = assay_name,
      target_class = class_name,
      vertex_order = row_number(),
      gate_source = source,
      n_control_points = nrow(class_points),
      n_trimmed_points = trimmed_n,
      min_hull_points = min_hull_points,
      polygon_expansion = polygon_expansion,
      .before = 1
    )
}

build_polygon_gates <- function(shared_droplets, geometry) {
  controls <- shared_droplets %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control"))

  bind_rows(lapply(mutation_order, function(assay_name) {
    bind_rows(lapply(class_names, function(class_name) {
      build_class_polygon(controls, geometry, assay_name, class_name)
    }))
  }))
}

close_polygon_paths <- function(polygons) {
  polygons %>%
    group_by(assay, target_class) %>%
    arrange(vertex_order, .by_group = TRUE) %>%
    group_modify(function(.x, .y) bind_rows(.x, slice(.x, 1))) %>%
    ungroup()
}

point_in_polygon <- function(x, y, polygon) {
  px <- polygon$ch1_amplitude
  py <- polygon$ch2_amplitude
  inside <- rep(FALSE, length(x))
  j <- length(px)

  for (i in seq_along(px)) {
    denominator <- py[j] - py[i]
    if (abs(denominator) < .Machine$double.eps) {
      denominator <- ifelse(denominator < 0, -.Machine$double.eps, .Machine$double.eps)
    }
    crosses <- ((py[i] > y) != (py[j] > y)) &
      (x < (px[j] - px[i]) * (y - py[i]) / denominator + px[i])
    crosses[is.na(crosses)] <- FALSE
    inside <- xor(inside, crosses)
    j <- i
  }

  inside
}

distance_matrix_for_geometry <- function(droplets, geometry, assay_name) {
  x <- as.matrix(droplets[, c("ch1_amplitude", "ch2_amplitude")])
  centres <- geometry$centres %>%
    filter(assay == assay_name, target_class %in% class_names) %>%
    arrange(match(target_class, class_names))

  distances <- sapply(seq_len(nrow(centres)), function(i) {
    class_name <- centres$target_class[[i]]
    mahalanobis(
      x,
      center = c(centres$centre_ch1[[i]], centres$centre_ch2[[i]]),
      cov = covariance_matrix_for(geometry$covariance_rows, assay_name, class_name)
    )
  })
  colnames(distances) <- centres$target_class
  distances
}

classify_polygon_droplets <- function(parsed, polygons, geometry) {
  droplets <- parsed$droplets
  assay_name <- parsed$assay
  assay_polygons <- polygons %>% filter(assay == assay_name)

  inside <- sapply(class_names, function(class_name) {
    polygon <- assay_polygons %>%
      filter(target_class == class_name) %>%
      arrange(vertex_order)
    point_in_polygon(droplets$ch1_amplitude, droplets$ch2_amplitude, polygon)
  })
  colnames(inside) <- class_names

  finite <- is.finite(droplets$ch1_amplitude) & is.finite(droplets$ch2_amplitude)
  inside[!finite, ] <- FALSE
  inside_count <- rowSums(inside)
  assigned <- rep("Rain", nrow(droplets))
  assigned[!finite] <- "Unclassified"

  one_match <- inside_count == 1L
  assigned[one_match] <- colnames(inside)[max.col(inside[one_match, , drop = FALSE])]

  multi_match <- inside_count > 1L
  if (any(multi_match)) {
    distances <- distance_matrix_for_geometry(droplets, geometry, assay_name)
    masked <- distances[multi_match, , drop = FALSE]
    masked[!inside[multi_match, , drop = FALSE]] <- Inf
    assigned[multi_match] <- colnames(masked)[max.col(-masked, ties.method = "first")]
  }

  droplets %>%
    mutate(
      assigned_target_class = factor(assigned, levels = target_class_levels),
      method_id = method_id,
      method_variant = method_variant,
      method_parameters_json = parameter_json
    )
}

failure_count_row <- function(parsed, message) {
  tibble(
    run_id = parsed$run_id,
    assay = parsed$assay,
    well = parsed$well,
    sample = parsed$sample,
    sample_type = parsed$sample_type,
    sample_key = parsed$sample_key,
    run_date = parsed$run_date,
    source_peak_data_json = parsed$peak_path,
    source_peak_metadata_json = parsed$metadata_path,
    method_id = method_id,
    method_variant = method_variant,
    n_total_droplets = nrow(parsed$droplets),
    n_accepted_droplets = NA_real_,
    n_nn = NA_real_,
    n_wt_only = NA_real_,
    n_mut_only = NA_real_,
    n_double_positive = NA_real_,
    n_rain = NA_real_,
    n_unclassified = NA_real_,
    n_mut_expected = NA_real_,
    n_wt_expected = NA_real_,
    classification_status = "failed",
    classification_message = message,
    method_parameters_json = parameter_json
  )
}

validation_rows <- function(classified) {
  classified %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
    count(
      method_id,
      method_variant,
      assay,
      sample_type,
      quanta_class = target_class,
      assigned_target_class,
      name = "droplets"
    ) %>%
    mutate(method_parameters_json = parameter_json)
}

diagnostic_sample_rows <- function(classified) {
  classified %>%
    group_by(assay, sample_type, assigned_target_class) %>%
    group_modify(function(.x, .y) {
      .x %>% slice_sample(n = min(nrow(.x), plot_sample_per_well))
    }) %>%
    ungroup() %>%
    select(
      method_id, method_variant, assay, run_id, well, sample, sample_type,
      ch1_amplitude, ch2_amplitude, target_class, assigned_target_class,
      method_parameters_json
    )
}

normalise_parameter_column <- function(rows) {
  lapply(rows, function(row) {
    row %>% mutate(method_parameters_json = as.character(method_parameters_json))
  })
}

save_polygon_gate_plots <- function(shared_droplets, polygons) {
  out_dir <- path_v2("plots", "individual", "polygon_gates")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  polygon_paths <- close_polygon_paths(polygons)
  plot_rows <- list()

  for (assay_name in mutation_order) {
    assay_controls <- shared_droplets %>%
      filter(assay == assay_name, sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
      group_by(sample_type, target_class) %>%
      group_modify(function(.x, .y) {
        .x %>% slice_sample(n = min(nrow(.x), 1000L))
      }) %>%
      ungroup()
    assay_polygons <- polygon_paths %>% filter(assay == assay_name)

    plot <- ggplot(assay_controls, aes(ch1_amplitude, ch2_amplitude, colour = target_class)) +
      geom_point(size = 0.16, alpha = 0.35, stroke = 0) +
      geom_path(
        data = assay_polygons,
        aes(ch1_amplitude, ch2_amplitude, group = target_class, colour = target_class),
        inherit.aes = FALSE,
        linewidth = 0.42
      ) +
      facet_wrap(~sample_type, nrow = 1) +
      coord_cartesian(
        xlim = c(0, quantile(assay_controls$ch1_amplitude, 0.999, na.rm = TRUE) * 1.05),
        ylim = c(0, quantile(assay_controls$ch2_amplitude, 0.999, na.rm = TRUE) * 1.05),
        expand = FALSE
      ) +
      labs(
        title = paste(assay_name, "control-derived polygon gates"),
        x = "Channel 1 amplitude",
        y = "Channel 2 amplitude",
        colour = "Class"
      ) +
      theme_bw(base_size = 8) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")

    svg_path <- file.path(out_dir, paste0(assay_name, "_polygon_gates.svg"))
    pdf_path <- file.path(out_dir, paste0(assay_name, "_polygon_gates.pdf"))
    grDevices::svg(svg_path, width = 8.5, height = 3.2, onefile = FALSE)
    print(plot)
    grDevices::dev.off()
    grDevices::pdf(pdf_path, width = 8.5, height = 3.2, useDingbats = FALSE)
    print(plot)
    grDevices::dev.off()

    plot_rows[[assay_name]] <- tibble(
      assay = assay_name,
      svg_path = svg_path,
      pdf_path = pdf_path
    )
  }

  bind_rows(plot_rows)
}

geometry <- readRDS(path_v2("models", "control_geometry", "control_geometry.rds"))
shared_droplets <- load_shared_droplets()
parsed_wells <- load_parsed_wells()
polygons <- build_polygon_gates(shared_droplets, geometry)

write_csv(
  tibble(
    method_id = method_id,
    method_variant = method_variant,
    min_hull_points = min_hull_points,
    polygon_expansion = polygon_expansion,
    ellipse_fallback_points = ellipse_fallback_points,
    random_seed = random_seed
  ),
  path_v2("tables", "polygon_parameter_grid.csv")
)
write_csv(polygons, path_v2("models", "polygon_gates", "polygon_vertices.csv"))

all_count_rows <- list()
all_validation_rows <- list()
all_sample_rows <- list()
all_status_rows <- list()

groups <- split(parsed_wells, vapply(parsed_wells, function(parsed) {
  paste(parsed$assay, parsed$run_id, sep = "::")
}, character(1)))

for (group_key in names(groups)) {
  group <- groups[[group_key]]
  assay_name <- group[[1]]$assay
  run_id <- group[[1]]$run_id
  elapsed <- system.time({
    status <- "ok"
    message <- NA_character_
    for (parsed in group) {
      result <- tryCatch(
        list(classified = classify_polygon_droplets(parsed, polygons, geometry), message = NA_character_),
        error = function(err) list(classified = NULL, message = conditionMessage(err))
      )

      if (is.null(result$classified)) {
        status <- "failed"
        message <- result$message
        all_count_rows <- append(all_count_rows, list(failure_count_row(parsed, result$message)))
        next
      }

      all_count_rows <- append(
        all_count_rows,
        list(droplet_count_row(
          parsed = parsed,
          method_id = method_id,
          method_variant = method_variant,
          target_class = result$classified$assigned_target_class,
          parameters = list(
            gate_source = "control_anchored_polygon",
            min_hull_points = min_hull_points,
            polygon_expansion = polygon_expansion
          )
        ))
      )
      all_validation_rows <- append(all_validation_rows, list(validation_rows(result$classified)))
      all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(result$classified)))
    }
  })

  all_status_rows <- append(all_status_rows, list(tibble(
    method_id = method_id,
    method_variant = method_variant,
    assay = assay_name,
    run_id = run_id,
    status = status,
    message = message,
    elapsed_seconds = max(0, unname(elapsed[["elapsed"]])),
    method_parameters_json = parameter_json
  )))
}

count_rows <- bind_rows(normalise_parameter_column(all_count_rows))
control_validation <- bind_rows(normalise_parameter_column(all_validation_rows))
plot_droplets <- bind_rows(all_sample_rows) %>%
  group_by(method_id, assay, sample_type, assigned_target_class) %>%
  group_modify(function(.x, .y) {
    .x %>% slice_sample(n = min(nrow(.x), plot_sample_per_group))
  }) %>%
  ungroup()
run_status <- bind_rows(all_status_rows)
plot_manifest <- save_polygon_gate_plots(shared_droplets, polygons)

write_csv(count_rows, path_v2("tables", "polygon_well_counts.csv"))
write_csv(control_validation, path_v2("tables", "polygon_control_validation.csv"))
write_csv(run_status, path_v2("tables", "polygon_run_status.csv"))
write_csv(plot_manifest, path_v2("tables", "polygon_plot_manifest.csv"))
saveRDS(plot_droplets, path_v2("data", "droplets", "polygon_control_hull_gates_plot_droplets.rds"))

e2e_checks <- tibble(
  check = c(
    "four_polygons_per_assay",
    "finite_polygon_vertices",
    "count_rows_for_complete_wells",
    "all_count_rows_ok",
    "control_validation_present",
    "plot_droplets_written",
    "polygon_plots_exist"
  ),
  passed = c(
    all(polygons %>% distinct(assay, target_class) %>% count(assay) %>% pull(n) == length(class_names)),
    all(is.finite(polygons$ch1_amplitude)) &&
      all(is.finite(polygons$ch2_amplitude)) &&
      all(polygons %>% count(assay, target_class) %>% pull(n) >= 3L),
    nrow(count_rows) == length(parsed_wells),
    all(count_rows$classification_status == "ok"),
    nrow(control_validation) > 0L,
    nrow(plot_droplets) > 0L &&
      file.exists(path_v2("data", "droplets", "polygon_control_hull_gates_plot_droplets.rds")),
    nrow(plot_manifest) == length(mutation_order) &&
      all(file.exists(c(plot_manifest$svg_path, plot_manifest$pdf_path))) &&
      all(file.info(c(plot_manifest$svg_path, plot_manifest$pdf_path))$size > 0)
  )
)
write_csv(e2e_checks, path_v2("tables", "polygon_e2e_checks.csv"))
stopifnot(all(e2e_checks$passed))

cat("polygon_vertices=", nrow(polygons), "\n", sep = "")
cat("well_count_rows=", nrow(count_rows), "\n", sep = "")
cat("control_validation_rows=", nrow(control_validation), "\n", sep = "")
cat("plot_droplet_rows=", nrow(plot_droplets), "\n", sep = "")
cat("e2e_checks_passed=", sum(e2e_checks$passed), "\n", sep = "")
