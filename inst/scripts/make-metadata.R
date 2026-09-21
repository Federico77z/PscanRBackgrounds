#!/usr/bin/env Rscript

if (!requireNamespace("ExperimentHubData", quietly = TRUE)) {
    stop(
        "Install ExperimentHubData before generating and validating metadata.",
        call. = FALSE
    )
}

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1]])
} else {
    "inst/scripts/make-metadata.R"
}
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
metadata <- data.frame(
    Title = "PscanR_backgrounds_v2",
    Description = paste(
        "Immutable ZIP archive containing 105 validated PscanR",
        "promoter-background score distributions for JASPAR 2020, 2022,",
        "and 2024, seven assemblies, and five promoter windows."
    ),
    BiocVersion = "3.24",
    Genome = NA_character_,
    SourceType = "Zip",
    SourceUrl = "https://doi.org/10.5281/zenodo.21821764",
    SourceVersion = "2",
    Species = NA_character_,
    TaxonomyId = NA_integer_,
    Coordinate_1_based = NA,
    DataProvider = "PscanR",
    Maintainer = "Federico Zambelli <federico.zambelli@unimi.it>",
    RDataClass = "character",
    DispatchClass = "FilePath",
    Location_Prefix = "https://zenodo.org/api/records/21821764/files/",
    RDataPath = "PscanR_backgrounds_v2.zip/content",
    Tags = "PscanR:GeneRegulation:MotifAnnotation:JASPAR:Promoter",
    stringsAsFactors = FALSE
)
utils::write.csv(
    metadata, file.path(root, "inst", "extdata", "metadata.csv"),
    row.names = FALSE, na = "NA"
)

ExperimentHubData::makeExperimentHubMetadata(root)
