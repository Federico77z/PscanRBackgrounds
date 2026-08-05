#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1]])
} else {
    "scripts/backgrounds.R"
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "background_pipeline.R"))
source(file.path(root, "R", "background_comparison.R"))

usage <- function(status = 0L) {
    cat(paste(
        "Usage: Rscript scripts/backgrounds.R MODE [OPTIONS]",
        "",
        "Modes:",
        "  audit      Validate existing artifacts and catalog entries",
        "  check      Compare current annotation-coordinate hashes with catalog",
        "  plan       Compute exact input hashes and write the regeneration plan",
        "  generate   Execute jobs from a new exact plan",
        "  validate   Validate all catalog artifacts and checksums",
        "  compare    Compare two complete background versions",
        "  calibrate  Run null calibration and promoter-universe comparison",
        "  all        Plan, generate, then validate",
        "",
        "Options:",
        "  --cores=N               BiocParallel workers (default: 1)",
        "  --repetitions=N         Calibration repetitions (default: 1000)",
        "  --minimum-set-size=N    Smallest gated calibration set (default: 50)",
        "  --assembly=A,B          Restrict assemblies",
        "  --jaspar=2020,2022      Restrict JASPAR releases",
        "  --window=450u_50d       Restrict promoter windows",
        "  --reference-version=N  Reference version for compare (default: 1)",
        "  --candidate-version=N  Candidate version for compare (default: 2)",
        "  --benchmark-dir=PATH    Existing foreground benchmark for compare",
        "  --output-dir=PATH       Comparison report directory",
        "  --force                 Regenerate selected jobs despite equal hashes",
        "  --help                  Show this help",
        sep = "\n"
    ))
    quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args == "--help")) usage(if (length(args)) 0L else 1L)
mode <- args[[1]]
args <- args[-1]
valid_modes <- c(
    "audit", "check", "plan", "generate", "validate", "compare",
    "calibrate", "all"
)
if (!mode %in% valid_modes) {
    cat("Unknown mode: ", mode, "\n", sep = "")
    usage(1L)
}

option_value <- function(name) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (!length(hit)) return(character())
    value <- sub(paste0("^--", name, "="), "", hit[[1]])
    strsplit(value, ",", fixed = TRUE)[[1]]
}

unknown <- args[!grepl(
    paste0(
        "^--(cores|repetitions|minimum-set-size|assembly|jaspar|window|",
        "reference-version|",
        "candidate-version|benchmark-dir|output-dir)=|^--force$"
    ), args
)]
if (length(unknown)) bg_stop("Unknown options: ", paste(unknown, collapse = ", "))
cores_value <- option_value("cores")
cores <- if (length(cores_value)) suppressWarnings(as.integer(cores_value[[1]])) else 1L
if (is.na(cores) || cores < 1L) bg_stop("--cores must be a positive integer")
repetitions_value <- option_value("repetitions")
repetitions <- if (length(repetitions_value)) {
    suppressWarnings(as.integer(repetitions_value[[1]]))
} else {
    1000L
}
if (is.na(repetitions) || repetitions < 10L) {
    bg_stop("--repetitions must be an integer of at least 10")
}
filters <- list(
    assembly = option_value("assembly"),
    jaspar = option_value("jaspar"),
    window = option_value("window")
)
force <- "--force" %in% args

positive_integer_option <- function(name, default) {
    value <- option_value(name)
    result <- if (length(value)) suppressWarnings(as.integer(value[[1]])) else default
    if (length(value) > 1L || is.na(result) || result < 1L) {
        bg_stop("--", name, " must be one positive integer")
    }
    result
}

path_option <- function(name, default = NULL, must_work = FALSE) {
    value <- option_value(name)
    if (!length(value)) return(default)
    if (length(value) != 1L || !nzchar(value[[1]])) {
        bg_stop("--", name, " must be one non-empty path")
    }
    path <- value[[1]]
    if (!grepl("^/", path)) path <- file.path(root, path)
    normalizePath(path, mustWork = must_work)
}

reference_version <- positive_integer_option("reference-version", 1L)
candidate_version <- positive_integer_option("candidate-version", 2L)
minimum_set_size <- positive_integer_option("minimum-set-size", 50L)
benchmark_dir <- path_option("benchmark-dir", must_work = TRUE)
comparison_output <- path_option(
    "output-dir",
    file.path(
        root, "reports",
        sprintf("comparison_v%d_v%d", reference_version, candidate_version)
    )
)

if (!nzchar(Sys.getenv("PSCANR_SOURCE", ""))) {
    sibling <- normalizePath(file.path(root, "..", "PscanR"), mustWork = FALSE)
    if (file.exists(file.path(sibling, "DESCRIPTION"))) {
        Sys.setenv(PSCANR_SOURCE = sibling)
    }
}

dir.create(file.path(root, "reports"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, ".cache"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "staging"), showWarnings = FALSE, recursive = TRUE)
config <- bg_read_config(root)
catalog <- bg_read_catalog(root)

write_plan <- function(plan) {
    path <- file.path(root, "reports", "generation_plan.tsv")
    utils::write.table(
        plan, path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    print(plan[, c(
        "action", "reason", "assembly", "upstream", "downstream",
        "jaspar_release", "background_version", "promoter_count", "motif_count"
    )], row.names = FALSE)
    bg_message(
        "Plan contains %d generation jobs and %d no-op jobs",
        sum(plan$action == "generate"), sum(plan$action == "skip")
    )
    plan
}

if (mode == "audit") {
    result <- bg_audit_existing(root, catalog)
    if (length(result$catalog) || any(nzchar(result$artifacts$problems))) {
        quit(save = "no", status = 1L)
    }
} else if (mode == "check") {
    bg_check_annotations(config, catalog, root, filters)
} else if (mode == "plan") {
    plan <- write_plan(bg_build_plan(config, catalog, root, filters, force))
} else if (mode == "generate") {
    plan <- write_plan(bg_build_plan(config, catalog, root, filters, force))
    catalog <- bg_generate(plan, catalog, root, cores)
} else if (mode == "validate") {
    problems <- bg_validate_catalog(catalog, root)
    if (length(problems)) {
        for (name in names(problems)) {
            cat(name, ": ", paste(problems[[name]], collapse = "; "), "\n", sep = "")
        }
        quit(save = "no", status = 1L)
    }
    bg_message("Validated %d catalog artifacts", nrow(catalog))
} else if (mode == "compare") {
    result <- bg_compare_versions(
        catalog = catalog,
        root = root,
        reference_version = reference_version,
        candidate_version = candidate_version,
        benchmark_dir = benchmark_dir,
        output_dir = comparison_output
    )
    print(result$checks, row.names = FALSE)
    bg_message("Comparison reports written to %s", comparison_output)
    if (any(!result$checks$passed)) quit(save = "no", status = 1L)
} else if (mode == "calibrate") {
    result <- bg_calibrate(
        config, root, filters, cores, repetitions,
        minimum_set_size = minimum_set_size
    )
    print(result$summary, row.names = FALSE)
    if (any(result$summary$assessment == "fail")) {
        quit(save = "no", status = 1L)
    }
} else if (mode == "all") {
    plan <- write_plan(bg_build_plan(config, catalog, root, filters, force))
    catalog <- bg_generate(plan, catalog, root, cores)
    problems <- bg_validate_catalog(catalog, root)
    if (length(problems)) bg_stop("Post-generation catalog validation failed")
    bg_message("Validated %d catalog artifacts", nrow(catalog))
}
