#!/usr/bin/env Rscript
# Assisted-by: OpenAI Codex and Claude Code (code, review and documentation).
# All changes were reviewed and tested by the authors.
# Build inst/extdata/metadata.csv: one ExperimentHub record per PscanR
# version-2 promoter background listed in inst/extdata/catalog_v2.tsv.
#
# The background files are hosted on Zenodo, one record per JASPAR release
# (35 files each), and each Hub record points at its file there. The same
# files are also archived, as a single ZIP, in the Zenodo record
# https://doi.org/10.5281/zenodo.21821764.
#
# Run from the package root:  Rscript inst/scripts/make-metadata.R

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1]])
} else {
    "inst/scripts/make-metadata.R"
}
root <- normalizePath(file.path(dirname(script_path), "..", ".."),
    mustWork = TRUE)

catalog <- utils::read.delim(
    file.path(root, "inst", "extdata", "catalog_v2.tsv"),
    colClasses = "character", na.strings = character()
)
catalog <- catalog[catalog$status == "validated", , drop = FALSE]
stopifnot(nrow(catalog) == 105L, all(catalog$background_version == "2"))

species <- c(
    hg38 = "Homo sapiens", hs1 = "Homo sapiens",
    mm10 = "Mus musculus", mm39 = "Mus musculus",
    dm6 = "Drosophila melanogaster",
    sacCer3 = "Saccharomyces cerevisiae",
    TAIR9 = "Arabidopsis thaliana"
)
taxonomy <- c(
    "Homo sapiens" = 9606L, "Mus musculus" = 10090L,
    "Drosophila melanogaster" = 7227L,
    "Saccharomyces cerevisiae" = 4932L,
    "Arabidopsis thaliana" = 3702L
)
annotation <- c(
    hg38 = "UCSC ncbiRefSeqCurated", mm10 = "UCSC ncbiRefSeqCurated",
    mm39 = "UCSC ncbiRefSeqCurated", dm6 = "UCSC ncbiRefSeqCurated",
    hs1 = "NCBI RefSeq GCF_009914755.1", sacCer3 = "UCSC sgdGene",
    TAIR9 = "TAIR9 GFF3"
)
provider <- c(
    hg38 = "JASPAR; UCSC", mm10 = "JASPAR; UCSC", mm39 = "JASPAR; UCSC",
    dm6 = "JASPAR; UCSC", hs1 = "JASPAR; NCBI", sacCer3 = "JASPAR; UCSC",
    TAIR9 = "JASPAR; TAIR"
)

# Zenodo records holding the individual version-2 files, by JASPAR release.
zenodo_record <- c("2020" = "23183695", "2022" = "23183713", "2024" = "23183720")
record <- unname(zenodo_record[catalog$jaspar_release])
stopifnot(!anyNA(record))

file <- basename(catalog$artifact)
stem <- sub("\\.psbg[0-9]+\\.txt$", "", file)
assembly <- catalog$assembly
organism <- unname(species[assembly])
stopifnot(!anyNA(organism))
window <- paste0(catalog$upstream, "u_", catalog$downstream, "d")

metadata <- data.frame(
    Title = paste0("PscanR_bg_v", catalog$background_version, "_", stem),
    Description = paste0(
        "PscanR version-2 promoter background for JASPAR ",
        catalog$jaspar_release, " CORE ", catalog$tax_group, " motifs (",
        catalog$motif_count, " matrices) scored on ",
        catalog$promoter_count, " unique ", assembly, " promoters from ",
        catalog$upstream, " bp upstream to ", catalog$downstream,
        " bp downstream of the TSS (", unname(annotation[assembly]),
        " annotation). Tab-separated text: motif identifier, number of ",
        "promoters, mean and standard deviation of the best PscanR score."
    ),
    BiocVersion = "3.24",
    Genome = assembly,
    SourceType = "TXT",
    SourceUrl = paste0("https://doi.org/10.5281/zenodo.", record),
    SourceVersion = catalog$background_version,
    Species = organism,
    TaxonomyId = unname(taxonomy[organism]),
    Coordinate_1_based = NA,
    DataProvider = unname(provider[assembly]),
    Maintainer = "Federico Zambelli <federico.zambelli@unimi.it>",
    RDataClass = "character",
    DispatchClass = "FilePath",
    RDataPath = paste0("records/", record, "/files/", file),
    Location_Prefix = "https://zenodo.org/",
    Tags = paste(
        "PscanR", "Promoter", "TFBS", paste0("JASPAR", catalog$jaspar_release),
        assembly, window, sep = ":"
    ),
    stringsAsFactors = FALSE
)
stopifnot(!anyDuplicated(metadata$Title), !anyDuplicated(metadata$RDataPath))

utils::write.csv(
    metadata, file.path(root, "inst", "extdata", "metadata.csv"),
    row.names = FALSE, na = "NA"
)

if (requireNamespace("ExperimentHubData", quietly = TRUE)) {
    invisible(ExperimentHubData::makeExperimentHubMetadata(root))
} else {
    message("ExperimentHubData is not installed; metadata not validated.")
}
