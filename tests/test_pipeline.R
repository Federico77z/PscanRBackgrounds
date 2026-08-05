script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1]])
} else {
    "tests/test_pipeline.R"
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "background_pipeline.R"))
source(file.path(root, "R", "background_comparison.R"))

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

prepared <- list(
    all_sequences = Biostrings::DNAStringSet(stats::setNames(
        c("AAAA", "CCCC", "GGGG"), c("tx1", "tx2", "tx3")
    )),
    snapshot = list(transcripts = data.frame(
        transcript_id = c("tx1", "tx2", "tx3"),
        gene_id = c("", "", ""), stringsAsFactors = FALSE
    ))
)
one_per_gene <- bg_one_promoter_per_gene(
    prepared, c(tx1 = "gene1", tx2 = "gene1", tx3 = "gene2")
)
stopifnot(length(one_per_gene) == 2L)

calibration_details <- data.frame(
    assembly = rep("test", 4L), set_size = rep(c(10L, 20L), each = 2L),
    ks_fdr = c(0.001, 0.5, 0.5, 0.5)
)
calibration_replicates <- data.frame(
    assembly = rep("test", 20L),
    set_size = rep(c(10L, 20L), each = 10L),
    fpr_005 = c(rep(0.10, 10L), rep(c(0.04, 0.06), 5L)),
    fpr_001 = c(rep(0.03, 10L), rep(c(0.005, 0.015), 5L))
)
calibration_summary <- bg_summarize_calibration(
    calibration_details, calibration_replicates, minimum_set_size = 20L
)
stopifnot(
    calibration_summary$assessment[calibration_summary$set_size == 10L] ==
        "small-set diagnostic",
    calibration_summary$assessment[calibration_summary$set_size == 20L] ==
        "pass"
)
failed_replicates <- calibration_replicates
failed_replicates$fpr_005[failed_replicates$set_size == 20L] <- 0.08
failed_summary <- bg_summarize_calibration(
    calibration_details, failed_replicates, minimum_set_size = 20L
)
stopifnot(
    failed_summary$assessment[failed_summary$set_size == 20L] == "fail"
)

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
catalog_keys <- bg_key(
    catalog$organism, catalog$assembly, catalog$upstream,
    catalog$downstream, catalog$jaspar_release
)
stopifnot(
    nrow(catalog) == 210L,
    length(unique(catalog_keys)) == 105L,
    all(table(catalog$background_version) == 105L),
    sum(catalog$latest) == 105L,
    all(catalog$background_version[catalog$latest] == 2L),
    all(catalog$status == "validated"),
    !anyDuplicated(catalog$artifact),
    !length(bg_validate_catalog(catalog, root))
)

key_info <- list(
    key = "synthetic", organism = "test", assembly = "test",
    upstream = 100L, downstream = 0L, jaspar_release = 2024L,
    reference_artifact = "reference.txt", reference_sha256 = "reference",
    candidate_artifact = "candidate.txt", candidate_sha256 = "candidate"
)
reference <- data.frame(
    motif_id = c("A", "B", "C"), background_size = 100L,
    mean = c(0.5, 0.6, 0.7), sd = c(0.05, 0.06, 0.07)
)
exact <- bg_compare_background_tables(reference, reference, key_info)
exact_checks <- bg_direct_acceptance(
    exact$summary, exact$detail, expected_pairs = 1L
)
stopifnot(all(exact_checks$passed))
plot_path <- tempfile(fileext = ".pdf")
bg_write_comparison_plots(
    exact$summary, exact$detail, benchmark = NULL, path = plot_path
)
stopifnot(file.info(plot_path)$size > 0L)
unlink(plot_path)

candidate <- reference
candidate$background_size <- 101L
candidate$mean <- candidate$mean + c(1e-5, -2e-5, 1e-5)
candidate$sd <- candidate$sd * c(1.001, 0.999, 1.0005)
changed <- bg_compare_background_tables(reference, candidate, key_info)
changed_checks <- bg_direct_acceptance(
    changed$summary, changed$detail, expected_pairs = 1L
)
stopifnot(all(changed_checks$passed))

mismatched <- reference[-1L, ]
mismatched$motif_id[[1]] <- "D"
motif_change <- bg_compare_background_tables(reference, mismatched, key_info)
motif_checks <- bg_direct_acceptance(
    motif_change$summary, motif_change$detail, expected_pairs = 1L
)
stopifnot(!motif_checks$passed[motif_checks$metric == "motif set differences"])

sentinel_reference <- data.frame(
    motif_id = "constant", background_size = 100L, mean = 1, sd = 1e-6
)
sentinel_candidate <- sentinel_reference
sentinel_candidate$sd <- 1e-5
sentinel <- bg_compare_background_tables(
    sentinel_reference, sentinel_candidate, key_info
)
sentinel_checks <- bg_direct_acceptance(
    sentinel$summary, sentinel$detail, expected_pairs = 1L
)
stopifnot(sentinel$detail$zero_variance_sentinel, all(sentinel_checks$passed))

benchmark_dir <- tempfile("background-benchmark-")
dir.create(file.path(benchmark_dir, "tables"), recursive = TRUE)
motif_count <- 60L
motif_id <- sprintf("M%03d", seq_len(motif_count))
n <- 20L
reference_bg <- data.frame(
    motif_id = motif_id, background_size = 100L,
    mean = 0.5 + seq_len(motif_count) / 1000,
    sd = rep(0.05, motif_count)
)
z <- seq(4, -2, length.out = motif_count)
foreground <- reference_bg$mean + z * reference_bg$sd / sqrt(n)
p_value <- stats::pnorm(z, lower.tail = FALSE)
benchmark_table <- data.frame(
    motif_id = motif_id, BG_AVG = reference_bg$mean,
    BG_STDEV = reference_bg$sd, FG_AVG = foreground, ZSCORE = z,
    P.VALUE = p_value, FDR = stats::p.adjust(p_value, "BH"),
    check.names = FALSE
)
utils::write.table(
    benchmark_table, file.path(benchmark_dir, "tables", "20.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
    data.frame(dataset = 20L, n_unique_sequences = n),
    file.path(benchmark_dir, "run_manifest.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
)
benchmark <- bg_compare_benchmark(
    benchmark_dir, reference_bg, reference_bg
)
stopifnot(
    nrow(benchmark$summary) == 1L,
    !nrow(benchmark$discordant),
    nzchar(benchmark$manifest_sha256),
    all(bg_benchmark_acceptance(benchmark)$passed)
)
unlink(benchmark_dir, recursive = TRUE)

cat("background pipeline tests passed\n")
