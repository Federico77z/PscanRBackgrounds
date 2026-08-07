#!/usr/bin/env Rscript

stop(paste(
    "Background generation requires the complete PscanR_backgrounds source",
    "repository. See scripts/backgrounds.R and README.md at",
    "https://github.com/Federico77z/PscanR_backgrounds. The published version-2",
    "archive is reproduced with: Rscript scripts/backgrounds.R zenodo",
    "--release-version=2"
), call. = FALSE)
