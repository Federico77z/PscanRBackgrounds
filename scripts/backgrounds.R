#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1]])
} else {
    "scripts/backgrounds.R"
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "background_pipeline.R"))

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
        "  calibrate  Run null calibration and promoter-universe comparison",
        "  all        Plan, generate, then validate",
        "",
        "Options:",
        "  --cores=N               BiocParallel workers (default: 1)",
        "  --repetitions=N         Calibration repetitions (default: 1000)",
        "  --assembly=A,B          Restrict assemblies",
        "  --jaspar=2020,2022      Restrict JASPAR releases",
        "  --window=450u_50d       Restrict promoter windows",
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
    "audit", "check", "plan", "generate", "validate", "calibrate", "all"
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
    "^--(cores|repetitions|assembly|jaspar|window)=|^--force$", args
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
} else if (mode == "calibrate") {
    result <- bg_calibrate(config, root, filters, cores, repetitions)
    print(result$summary, row.names = FALSE)
    if (any(!result$summary$passed)) quit(save = "no", status = 1L)
} else if (mode == "all") {
    plan <- write_plan(bg_build_plan(config, catalog, root, filters, force))
    catalog <- bg_generate(plan, catalog, root, cores)
    problems <- bg_validate_catalog(catalog, root)
    if (length(problems)) bg_stop("Post-generation catalog validation failed")
    bg_message("Validated %d catalog artifacts", nrow(catalog))
}
