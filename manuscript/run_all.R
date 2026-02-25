# Entry point for manuscript figure/table generation.
# Add calls to figure/table scripts as they are implemented.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  stop("Could not determine script path. Run with: Rscript manuscript/run_all.R")
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

message("Project root: ", project_root)
message("No manuscript scripts are wired yet. Add source(...) calls in manuscript/run_all.R.")

# Example:
# source(file.path(project_root, "manuscript", "scripts", "fig_01_template.R"))
# source(file.path(project_root, "manuscript", "scripts", "table_01_template.R"))
