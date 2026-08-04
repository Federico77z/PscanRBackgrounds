script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1]])
} else {
    "tests/test_pipeline.R"
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "background_pipeline.R"))

config <- bg_read_config(root)
stopifnot(
    nrow(config$organisms) == 7L,
    nrow(config$windows) == 5L,
    nrow(config$jaspar) == 3L,
    nrow(config$organisms) * nrow(config$windows) * nrow(config$jaspar) == 105L
)

stopifnot(identical(
    bg_filename(2024L, "hg38", 950L, 50L, "ucsc", 2L),
    "J2024_hg38_950u_50d_UCSC.psbg2.txt"
))
stopifnot(identical(
    bg_filename(2024L, "TAIR9", 1000L, 0L, "gff", 2L),
    "J2024_TAIR9_1000u_0d_TAIR.psbg2.txt"
))

transcripts <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4L),
    ranges = IRanges::IRanges(seq_len(4L), width = 1L)
)
S4Vectors::mcols(transcripts)$tx_id <- seq_len(4L)
S4Vectors::mcols(transcripts)$tx_name <- c(
    "NM_000001.1", "XM_000002.1", NA_character_, ""
)
curated <- bg_name_transcripts(transcripts, curated_refseq = TRUE)
stopifnot(identical(names(curated), "NM_000001.1"))
all_transcripts <- bg_name_transcripts(transcripts)
stopifnot(identical(
    names(all_transcripts),
    c("NM_000001.1", "XM_000002.1", "tx-3", "tx-4")
))

temporary <- tempfile(fileext = ".txt")
writeLines(c(
    "[SHORT TFBS MATRIX]",
    "MA0001.1\t10\t0.5\t0.1",
    "MA0002.1\t10\t0.6\t0.2"
), temporary)
stopifnot(length(bg_validate_artifact(temporary, 2L)) == 0L)
stopifnot(length(bg_validate_artifact(temporary, 3L)) == 1L)
unlink(temporary)

catalog <- bg_read_catalog(root)
stopifnot(
    nrow(catalog) == 105L,
    sum(catalog$latest) == 105L,
    all(catalog$status == "validated"),
    !anyDuplicated(catalog$artifact),
    !length(bg_validate_catalog(catalog, root))
)

cat("background pipeline tests passed\n")
