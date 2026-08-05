BG_COMPARISON_LIMITS <- list(
    stable_mean_tolerance = 1e-12,
    stable_sd_tolerance = 1e-12,
    correlation_minimum = 0.9999,
    normalized_mean_shift_maximum = 0.01,
    sd_ratio_minimum = 0.95,
    sd_ratio_maximum = 1.05,
    sd_ratio_p01_minimum = 0.98,
    sd_ratio_p99_maximum = 1.02,
    projected_z_n1000_maximum = 0.25,
    rank_correlation_minimum = 0.999,
    top10_overlap_minimum = 0.80,
    top25_overlap_minimum = 0.90,
    top50_overlap_minimum = 0.90,
    benchmark_z_delta_maximum = 0.10,
    fdr_discordance_maximum = 0.005,
    fdr_boundary_distance_maximum = 0.01
)

bg_read_background_table <- function(path) {
    table <- utils::read.delim(
        path, skip = 1L, header = FALSE, stringsAsFactors = FALSE
    )
    if (ncol(table) != 4L) bg_stop("Expected four columns in ", path)
    names(table) <- c("motif_id", "background_size", "mean", "sd")
    table
}

bg_safe_correlation <- function(x, y) {
    if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) {
        return(NA_real_)
    }
    stats::cor(x, y)
}

bg_compare_background_tables <- function(reference, candidate, key_info) {
    names(reference) <- c("motif_id", "reference_size", "reference_mean",
                          "reference_sd")
    names(candidate) <- c("motif_id", "candidate_size", "candidate_mean",
                          "candidate_sd")
    detail <- merge(
        reference, candidate, by = "motif_id", all = TRUE, sort = TRUE
    )
    detail$motif_status <- ifelse(
        is.na(detail$reference_size), "added",
        ifelse(is.na(detail$candidate_size), "removed", "shared")
    )
    shared <- detail$motif_status == "shared"
    detail$mean_delta <- detail$candidate_mean - detail$reference_mean
    detail$sd_delta <- detail$candidate_sd - detail$reference_sd
    detail$normalized_mean_shift <- detail$mean_delta / detail$candidate_sd
    detail$sd_ratio <- detail$reference_sd / detail$candidate_sd
    for (size in c(10L, 20L, 50L, 100L, 200L, 1000L)) {
        name <- paste0("projected_null_z_n", size)
        detail[[name]] <- -detail$normalized_mean_shift * sqrt(size)
    }
    detail$zero_variance_sentinel <- shared &
        abs(detail$reference_mean - 1) <= 1e-12 &
        abs(detail$candidate_mean - 1) <= 1e-12 &
        abs(detail$reference_sd - 1e-6) <= 1e-12 &
        abs(detail$candidate_sd - 1e-5) <= 1e-12
    for (name in names(key_info)) detail[[name]] <- key_info[[name]]

    current <- detail[shared, , drop = FALSE]
    summary <- data.frame(
        key = key_info$key,
        organism = key_info$organism,
        assembly = key_info$assembly,
        upstream = key_info$upstream,
        downstream = key_info$downstream,
        jaspar_release = key_info$jaspar_release,
        reference_artifact = key_info$reference_artifact,
        reference_sha256 = key_info$reference_sha256,
        candidate_artifact = key_info$candidate_artifact,
        candidate_sha256 = key_info$candidate_sha256,
        reference_size = unique(reference$reference_size),
        candidate_size = unique(candidate$candidate_size),
        reference_motifs = nrow(reference),
        candidate_motifs = nrow(candidate),
        shared_motifs = sum(shared),
        added_motifs = sum(detail$motif_status == "added"),
        removed_motifs = sum(detail$motif_status == "removed"),
        mean_correlation = bg_safe_correlation(
            current$reference_mean, current$candidate_mean
        ),
        sd_correlation = bg_safe_correlation(
            current$reference_sd, current$candidate_sd
        ),
        mean_absolute_delta = mean(abs(current$mean_delta)),
        mean_maximum_delta = max(abs(current$mean_delta)),
        sd_absolute_delta = mean(abs(current$sd_delta)),
        sd_maximum_delta = max(abs(current$sd_delta)),
        normalized_mean_shift_p95 = unname(stats::quantile(
            abs(current$normalized_mean_shift), 0.95
        )),
        normalized_mean_shift_maximum = max(
            abs(current$normalized_mean_shift)
        ),
        sd_ratio_p01 = unname(stats::quantile(current$sd_ratio, 0.01)),
        sd_ratio_median = stats::median(current$sd_ratio),
        sd_ratio_p99 = unname(stats::quantile(current$sd_ratio, 0.99)),
        sd_ratio_minimum = min(current$sd_ratio),
        sd_ratio_maximum = max(current$sd_ratio),
        zero_variance_sentinels = sum(current$zero_variance_sentinel),
        stringsAsFactors = FALSE
    )
    summary$promoter_count_delta <-
        summary$candidate_size - summary$reference_size
    summary$promoter_count_relative_delta <-
        summary$promoter_count_delta / summary$reference_size
    list(summary = summary, detail = detail)
}

bg_check_row <- function(scope, metric, observed, criterion, passed) {
    data.frame(
        scope = scope, metric = metric, observed = as.character(observed),
        criterion = criterion, passed = isTRUE(passed),
        stringsAsFactors = FALSE
    )
}

bg_max_or_zero <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else 0
}

bg_min_or_one <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else 1
}

bg_direct_acceptance <- function(key_summary, motif_detail,
                                 expected_pairs = 105L,
                                 limits = BG_COMPARISON_LIMITS) {
    stable_keys <- key_summary$reference_size == key_summary$candidate_size
    stable <- motif_detail[
        motif_detail$key %in% key_summary$key[stable_keys] &
            motif_detail$motif_status == "shared",
        , drop = FALSE
    ]
    stable_regular <- stable[!stable$zero_variance_sentinel, , drop = FALSE]
    changed_keys <- key_summary$reference_size != key_summary$candidate_size
    changed <- motif_detail[
        motif_detail$key %in% key_summary$key[changed_keys] &
            motif_detail$motif_status == "shared",
        , drop = FALSE
    ]
    changed_correlations <- c(
        key_summary$mean_correlation[changed_keys],
        key_summary$sd_correlation[changed_keys]
    )
    changed_correlations <- changed_correlations[
        is.finite(changed_correlations)
    ]
    ratio <- changed$sd_ratio[is.finite(changed$sd_ratio)]
    ratio_p01 <- if (length(ratio)) {
        unname(stats::quantile(ratio, 0.01))
    } else 1
    ratio_p99 <- if (length(ratio)) {
        unname(stats::quantile(ratio, 0.99))
    } else 1

    do.call(rbind, list(
        bg_check_row(
            "profiles", "paired background keys", nrow(key_summary),
            paste0("= ", expected_pairs), nrow(key_summary) == expected_pairs
        ),
        bg_check_row(
            "profiles", "motif set differences",
            sum(key_summary$added_motifs + key_summary$removed_motifs),
            "= 0",
            all(key_summary$added_motifs + key_summary$removed_motifs == 0)
        ),
        bg_check_row(
            "stable profiles", "maximum mean delta",
            bg_max_or_zero(abs(stable_regular$mean_delta)),
            paste0("<= ", limits$stable_mean_tolerance),
            bg_max_or_zero(abs(stable_regular$mean_delta)) <=
                limits$stable_mean_tolerance
        ),
        bg_check_row(
            "stable profiles", "maximum SD delta excluding sentinels",
            bg_max_or_zero(abs(stable_regular$sd_delta)),
            paste0("<= ", limits$stable_sd_tolerance),
            bg_max_or_zero(abs(stable_regular$sd_delta)) <=
                limits$stable_sd_tolerance
        ),
        bg_check_row(
            "changed profiles", "minimum mean/SD correlation",
            bg_min_or_one(changed_correlations),
            paste0(">= ", limits$correlation_minimum),
            bg_min_or_one(changed_correlations) >= limits$correlation_minimum
        ),
        bg_check_row(
            "changed profiles", "maximum normalized mean shift",
            bg_max_or_zero(abs(changed$normalized_mean_shift)),
            paste0("<= ", limits$normalized_mean_shift_maximum),
            bg_max_or_zero(abs(changed$normalized_mean_shift)) <=
                limits$normalized_mean_shift_maximum
        ),
        bg_check_row(
            "changed profiles", "minimum SD ratio",
            bg_min_or_one(ratio), paste0(">= ", limits$sd_ratio_minimum),
            bg_min_or_one(ratio) >= limits$sd_ratio_minimum
        ),
        bg_check_row(
            "changed profiles", "maximum SD ratio",
            bg_max_or_zero(ratio), paste0("<= ", limits$sd_ratio_maximum),
            bg_max_or_zero(ratio) <= limits$sd_ratio_maximum
        ),
        bg_check_row(
            "changed profiles", "SD ratio first percentile", ratio_p01,
            paste0(">= ", limits$sd_ratio_p01_minimum),
            ratio_p01 >= limits$sd_ratio_p01_minimum
        ),
        bg_check_row(
            "changed profiles", "SD ratio 99th percentile", ratio_p99,
            paste0("<= ", limits$sd_ratio_p99_maximum),
            ratio_p99 <= limits$sd_ratio_p99_maximum
        ),
        bg_check_row(
            "changed profiles", "maximum projected null Z shift at n=1000",
            bg_max_or_zero(abs(changed$projected_null_z_n1000)),
            paste0("<= ", limits$projected_z_n1000_maximum),
            bg_max_or_zero(abs(changed$projected_null_z_n1000)) <=
                limits$projected_z_n1000_maximum
        )
    ))
}

bg_top_overlap <- function(reference_p, candidate_p, count) {
    count <- min(count, length(reference_p))
    reference <- order(reference_p)[seq_len(count)]
    candidate <- order(candidate_p)[seq_len(count)]
    length(intersect(reference, candidate)) / count
}

bg_compare_benchmark <- function(benchmark_dir, reference_bg, candidate_bg) {
    manifest_path <- file.path(benchmark_dir, "run_manifest.tsv")
    if (!file.exists(manifest_path)) {
        bg_stop("Missing benchmark manifest: ", manifest_path)
    }
    manifest <- utils::read.delim(manifest_path, stringsAsFactors = FALSE)
    required_manifest <- c("dataset", "n_unique_sequences")
    if (!all(required_manifest %in% names(manifest))) {
        bg_stop("Benchmark manifest lacks dataset or n_unique_sequences")
    }
    if (anyDuplicated(manifest$dataset) ||
            any(!is.finite(manifest$n_unique_sequences)) ||
            any(manifest$n_unique_sequences < 1L)) {
        bg_stop("Benchmark manifest has invalid or duplicate datasets")
    }
    summaries <- list()
    discordant <- list()
    discordant_index <- 0L
    for (i in seq_len(nrow(manifest))) {
        table_path <- file.path(
            benchmark_dir, "tables", paste0(manifest$dataset[[i]], ".tsv")
        )
        table <- utils::read.delim(
            table_path, stringsAsFactors = FALSE, check.names = FALSE
        )
        required <- c(
            "motif_id", "BG_AVG", "BG_STDEV", "FG_AVG", "ZSCORE",
            "P.VALUE", "FDR"
        )
        if (!all(required %in% names(table))) {
            bg_stop("Invalid benchmark result table: ", table_path)
        }
        joined <- merge(
            table, reference_bg, by = "motif_id", sort = FALSE,
            suffixes = c("", "_reference")
        )
        joined <- merge(
            joined, candidate_bg, by = "motif_id", sort = FALSE,
            suffixes = c("_reference", "_candidate")
        )
        if (nrow(joined) != nrow(table) ||
                !setequal(joined$motif_id, table$motif_id)) {
            bg_stop(
                "Benchmark motifs do not match both backgrounds: ", table_path
            )
        }
        n <- manifest$n_unique_sequences[[i]]
        reference_z <- (joined$FG_AVG - joined$mean_reference) /
            (joined$sd_reference / sqrt(n))
        candidate_z <- (joined$FG_AVG - joined$mean_candidate) /
            (joined$sd_candidate / sqrt(n))
        candidate_p <- stats::pnorm(candidate_z, lower.tail = FALSE)
        candidate_fdr <- stats::p.adjust(candidate_p, method = "BH")
        reference_rank <- rank(joined$P.VALUE, ties.method = "min")
        candidate_rank <- rank(candidate_p, ties.method = "min")
        changed <- (joined$FDR < 0.05) != (candidate_fdr < 0.05)
        if (any(changed)) {
            discordant_index <- discordant_index + 1L
            discordant[[discordant_index]] <- data.frame(
                dataset = manifest$dataset[[i]], n = n,
                motif_id = joined$motif_id[changed],
                reference_z = joined$ZSCORE[changed],
                candidate_z = candidate_z[changed],
                reference_fdr = joined$FDR[changed],
                candidate_fdr = candidate_fdr[changed],
                direction = ifelse(
                    candidate_fdr[changed] < 0.05,
                    "candidate_only", "reference_only"
                ),
                stringsAsFactors = FALSE
            )
        }
        summaries[[i]] <- data.frame(
            dataset = manifest$dataset[[i]], n = n,
            table_sha256 = bg_hash_file(table_path),
            motifs = nrow(joined),
            reference_background_maximum_error = max(abs(c(
                joined$BG_AVG - joined$mean_reference,
                joined$BG_STDEV - joined$sd_reference
            ))),
            reference_recalculation_maximum_error = max(
                abs(reference_z - joined$ZSCORE)
            ),
            maximum_absolute_z_delta = max(
                abs(candidate_z - joined$ZSCORE)
            ),
            z_delta_p95 = unname(stats::quantile(
                abs(candidate_z - joined$ZSCORE), 0.95
            )),
            rank_spearman = stats::cor(
                reference_rank, candidate_rank, method = "spearman"
            ),
            top10_overlap = bg_top_overlap(
                joined$P.VALUE, candidate_p, 10L
            ),
            top25_overlap = bg_top_overlap(
                joined$P.VALUE, candidate_p, 25L
            ),
            top50_overlap = bg_top_overlap(
                joined$P.VALUE, candidate_p, 50L
            ),
            reference_significant = sum(joined$FDR < 0.05),
            candidate_significant = sum(candidate_fdr < 0.05),
            discordant_calls = sum(changed),
            maximum_absolute_fdr_delta = max(
                abs(candidate_fdr - joined$FDR)
            ),
            stringsAsFactors = FALSE
        )
    }
    summary <- do.call(rbind, summaries)
    discordant <- if (length(discordant)) {
        do.call(rbind, discordant)
    } else {
        data.frame(
            dataset = integer(), n = integer(), motif_id = character(),
            reference_z = numeric(), candidate_z = numeric(),
            reference_fdr = numeric(), candidate_fdr = numeric(),
            direction = character(), stringsAsFactors = FALSE
        )
    }
    list(
        summary = summary, discordant = discordant,
        manifest_sha256 = bg_hash_file(manifest_path)
    )
}

bg_benchmark_acceptance <- function(result,
                                    limits = BG_COMPARISON_LIMITS) {
    summary <- result$summary
    discordant <- result$discordant
    total <- sum(summary$motifs)
    discordance_rate <- nrow(discordant) / total
    boundary_distance <- if (nrow(discordant)) {
        max(abs(c(
            discordant$reference_fdr, discordant$candidate_fdr
        ) - 0.05))
    } else 0
    do.call(rbind, list(
        bg_check_row(
            "benchmark", "reference background value error",
            max(summary$reference_background_maximum_error), "<= 1e-12",
            max(summary$reference_background_maximum_error) <= 1e-12
        ),
        bg_check_row(
            "benchmark", "reference Z recalculation error",
            max(summary$reference_recalculation_maximum_error), "<= 1e-10",
            max(summary$reference_recalculation_maximum_error) <= 1e-10
        ),
        bg_check_row(
            "benchmark", "minimum rank Spearman",
            min(summary$rank_spearman),
            paste0(">= ", limits$rank_correlation_minimum),
            min(summary$rank_spearman) >= limits$rank_correlation_minimum
        ),
        bg_check_row(
            "benchmark", "minimum top-10 overlap",
            min(summary$top10_overlap),
            paste0(">= ", limits$top10_overlap_minimum),
            min(summary$top10_overlap) >= limits$top10_overlap_minimum
        ),
        bg_check_row(
            "benchmark", "minimum top-25 overlap",
            min(summary$top25_overlap),
            paste0(">= ", limits$top25_overlap_minimum),
            min(summary$top25_overlap) >= limits$top25_overlap_minimum
        ),
        bg_check_row(
            "benchmark", "minimum top-50 overlap",
            min(summary$top50_overlap),
            paste0(">= ", limits$top50_overlap_minimum),
            min(summary$top50_overlap) >= limits$top50_overlap_minimum
        ),
        bg_check_row(
            "benchmark", "maximum absolute Z delta",
            max(summary$maximum_absolute_z_delta),
            paste0("<= ", limits$benchmark_z_delta_maximum),
            max(summary$maximum_absolute_z_delta) <=
                limits$benchmark_z_delta_maximum
        ),
        bg_check_row(
            "benchmark", "FDR classification discordance rate",
            discordance_rate,
            paste0("<= ", limits$fdr_discordance_maximum),
            discordance_rate <= limits$fdr_discordance_maximum
        ),
        bg_check_row(
            "benchmark", "maximum discordant FDR boundary distance",
            boundary_distance,
            paste0("<= ", limits$fdr_boundary_distance_maximum),
            boundary_distance <= limits$fdr_boundary_distance_maximum
        )
    ))
}

bg_comparison_catalog_pairs <- function(catalog, reference_version,
                                        candidate_version) {
    reference <- catalog[
        catalog$status == "validated" &
            catalog$background_version == reference_version,
        , drop = FALSE
    ]
    candidate <- catalog[
        catalog$status == "validated" &
            catalog$background_version == candidate_version,
        , drop = FALSE
    ]
    reference$key <- bg_key(
        reference$organism, reference$assembly, reference$upstream,
        reference$downstream, reference$jaspar_release
    )
    candidate$key <- bg_key(
        candidate$organism, candidate$assembly, candidate$upstream,
        candidate$downstream, candidate$jaspar_release
    )
    if (anyDuplicated(reference$key) || anyDuplicated(candidate$key)) {
        bg_stop("A comparison version contains duplicate background keys")
    }
    if (!setequal(reference$key, candidate$key)) {
        missing_candidate <- setdiff(reference$key, candidate$key)
        missing_reference <- setdiff(candidate$key, reference$key)
        bg_stop(
            "Background version keys do not match. Missing candidate: ",
            paste(missing_candidate, collapse = ", "),
            "; missing reference: ",
            paste(missing_reference, collapse = ", ")
        )
    }
    candidate <- candidate[match(reference$key, candidate$key), ]
    list(reference = reference, candidate = candidate)
}

bg_write_comparison_plots <- function(key_summary, motif_detail, benchmark,
                                      path) {
    grDevices::pdf(path, width = 10, height = 7, onefile = TRUE)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit({
        graphics::par(old_par)
        grDevices::dev.off()
    }, add = TRUE)

    colors <- as.integer(factor(key_summary$assembly))
    graphics::plot(
        key_summary$reference_size, key_summary$candidate_size,
        col = colors, pch = 19, xlab = "Version 1 promoter count",
        ylab = "Version 2 promoter count", main = "Promoter universe sizes"
    )
    graphics::abline(0, 1, lty = 2)
    graphics::legend(
        "topleft", legend = levels(factor(key_summary$assembly)),
        col = seq_along(levels(factor(key_summary$assembly))), pch = 19,
        cex = 0.8
    )

    shared <- motif_detail$motif_status == "shared"
    graphics::plot(
        motif_detail$reference_mean[shared],
        motif_detail$candidate_mean[shared], pch = ".",
        xlab = "Version 1 mean", ylab = "Version 2 mean",
        main = "Background means"
    )
    graphics::abline(0, 1, col = "red", lty = 2)

    graphics::plot(
        motif_detail$reference_sd[shared], motif_detail$candidate_sd[shared],
        log = "xy", pch = ".", xlab = "Version 1 SD",
        ylab = "Version 2 SD", main = "Background standard deviations"
    )
    graphics::abline(0, 1, col = "red", lty = 2)

    changed_keys <- key_summary$key[
        key_summary$reference_size != key_summary$candidate_size
    ]
    changed <- shared & motif_detail$key %in% changed_keys
    if (any(changed)) {
        graphics::hist(
            motif_detail$normalized_mean_shift[changed], breaks = 60,
            xlab = "(version 2 mean - version 1 mean) / version 2 SD",
            main = "Normalized mean shifts for changed promoter universes"
        )
    } else {
        graphics::plot.new()
        graphics::title("Normalized mean shifts")
        graphics::text(0.5, 0.5, "No promoter-universe changes")
    }

    if (!is.null(benchmark)) {
        graphics::plot(
            benchmark$summary$n, benchmark$summary$maximum_absolute_z_delta,
            pch = 19, xlab = "Unique foreground promoters",
            ylab = "Maximum absolute Z-score delta",
            main = "Downstream benchmark Z-score impact"
        )
        graphics::abline(
            h = BG_COMPARISON_LIMITS$benchmark_z_delta_maximum,
            col = "red", lty = 2
        )
        graphics::plot(
            benchmark$summary$n, benchmark$summary$rank_spearman,
            pch = 19, xlab = "Unique foreground promoters",
            ylab = "Spearman correlation", ylim = c(0.999, 1),
            main = "Downstream motif-rank concordance"
        )
        graphics::abline(
            h = BG_COMPARISON_LIMITS$rank_correlation_minimum,
            col = "red", lty = 2
        )
    }
}

bg_format_number <- function(x, digits = 6L) {
    format(signif(x, digits), scientific = TRUE, trim = TRUE)
}

bg_write_comparison_summary <- function(result, reference_version,
                                        candidate_version, path) {
    key_summary <- result$key_summary
    motif_detail <- result$motif_detail
    changed <- key_summary$reference_size != key_summary$candidate_size
    shared_changed <- motif_detail$motif_status == "shared" &
        motif_detail$key %in% key_summary$key[changed]
    lines <- c(
        sprintf(
            "# Background version %d versus version %d",
            reference_version, candidate_version
        ),
        "",
        sprintf(
            "Overall result: **%s**.",
            if (all(result$checks$passed)) "PASS" else "FAIL"
        ),
        "",
        "## Profile comparison",
        "",
        sprintf("- Paired background keys: %d.", nrow(key_summary)),
        sprintf(
            "- Keys with changed promoter counts: %d; unchanged: %d.",
            sum(changed), sum(!changed)
        ),
        sprintf(
            "- Motif additions/removals: %d/%d.",
            sum(key_summary$added_motifs), sum(key_summary$removed_motifs)
        ),
        sprintf(
            "- Maximum normalized mean shift in changed universes: %s SD.",
            bg_format_number(bg_max_or_zero(abs(
                motif_detail$normalized_mean_shift[shared_changed]
            )))
        ),
        sprintf(
            "- Maximum projected null Z-score shift at n=1000: %s.",
            bg_format_number(bg_max_or_zero(abs(
                motif_detail$projected_null_z_n1000[shared_changed]
            )))
        ),
        sprintf(
            "- Documented zero-variance SD sentinel changes: %d motifs.",
            sum(motif_detail$zero_variance_sentinel, na.rm = TRUE)
        )
    )
    if (!is.null(result$benchmark)) {
        benchmark <- result$benchmark
        total <- sum(benchmark$summary$motifs)
        lines <- c(
            lines, "", "## Downstream benchmark", "",
            sprintf(
                "- Independent foreground datasets: %d.",
                nrow(benchmark$summary)
            ),
            sprintf(
                "- Benchmark manifest SHA-256: `%s`.",
                benchmark$manifest_sha256
            ),
            sprintf(
                "- Maximum absolute Z-score delta: %s.",
                bg_format_number(max(
                    benchmark$summary$maximum_absolute_z_delta
                ))
            ),
            sprintf(
                "- Minimum rank Spearman correlation: %s.",
                bg_format_number(min(benchmark$summary$rank_spearman))
            ),
            sprintf(
                "- FDR classification changes: %d of %d motif/dataset pairs (%s).",
                nrow(benchmark$discordant), total,
                bg_format_number(nrow(benchmark$discordant) / total)
            )
        )
    }
    lines <- c(
        lines, "", "## Acceptance checks", "",
        "| Scope | Metric | Observed | Criterion | Passed |",
        "|---|---|---:|---:|:---:|",
        apply(result$checks, 1L, function(row) {
            sprintf(
                "| %s | %s | %s | %s | %s |",
                row[["scope"]], row[["metric"]], row[["observed"]],
                row[["criterion"]], if (row[["passed"]] == "TRUE") "yes" else "no"
            )
        })
    )
    writeLines(lines, path)
}

bg_write_comparison_reports <- function(result, reference_version,
                                        candidate_version, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
        result$key_summary, file.path(output_dir, "key_summary.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    connection <- gzfile(
        file.path(output_dir, "motif_deltas.tsv.gz"), open = "wt"
    )
    on.exit(close(connection), add = TRUE)
    utils::write.table(
        result$motif_detail, connection, sep = "\t", quote = FALSE,
        row.names = FALSE
    )
    close(connection)
    on.exit(NULL, add = FALSE)
    utils::write.table(
        result$checks, file.path(output_dir, "acceptance_checks.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    if (!is.null(result$benchmark)) {
        utils::write.table(
            result$benchmark$summary,
            file.path(output_dir, "benchmark_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        utils::write.table(
            result$benchmark$discordant,
            file.path(output_dir, "benchmark_fdr_discordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
    bg_write_comparison_plots(
        result$key_summary, result$motif_detail, result$benchmark,
        file.path(output_dir, "comparison_plots.pdf")
    )
    bg_write_comparison_summary(
        result, reference_version, candidate_version,
        file.path(output_dir, "SUMMARY.md")
    )
    invisible(output_dir)
}

bg_compare_versions <- function(catalog, root, reference_version = 1L,
                                candidate_version = 2L,
                                benchmark_dir = NULL,
                                output_dir = file.path(
                                    root, "reports", "comparison_v1_v2"
                                )) {
    pairs <- bg_comparison_catalog_pairs(
        catalog, reference_version, candidate_version
    )
    selected_artifacts <- rbind(pairs$reference, pairs$candidate)
    invalid_hash <- vapply(seq_len(nrow(selected_artifacts)), function(i) {
        row <- selected_artifacts[i, , drop = FALSE]
        path <- file.path(root, row$artifact[[1]])
        !file.exists(path) ||
            !identical(bg_hash_file(path), row$artifact_sha256[[1]])
    }, logical(1))
    if (any(invalid_hash)) {
        bg_stop(
            "Missing or checksum-invalid comparison artifacts: ",
            paste(selected_artifacts$artifact[invalid_hash], collapse = ", ")
        )
    }
    comparisons <- lapply(seq_len(nrow(pairs$reference)), function(i) {
        reference_row <- pairs$reference[i, , drop = FALSE]
        candidate_row <- pairs$candidate[i, , drop = FALSE]
        key_info <- list(
            key = reference_row$key[[1]],
            organism = reference_row$organism[[1]],
            assembly = reference_row$assembly[[1]],
            upstream = reference_row$upstream[[1]],
            downstream = reference_row$downstream[[1]],
            jaspar_release = reference_row$jaspar_release[[1]],
            reference_artifact = reference_row$artifact[[1]],
            reference_sha256 = reference_row$artifact_sha256[[1]],
            candidate_artifact = candidate_row$artifact[[1]],
            candidate_sha256 = candidate_row$artifact_sha256[[1]]
        )
        bg_compare_background_tables(
            bg_read_background_table(file.path(root, reference_row$artifact)),
            bg_read_background_table(file.path(root, candidate_row$artifact)),
            key_info
        )
    })
    key_summary <- do.call(rbind, lapply(comparisons, `[[`, "summary"))
    motif_detail <- do.call(rbind, lapply(comparisons, `[[`, "detail"))
    row.names(key_summary) <- NULL
    row.names(motif_detail) <- NULL
    checks <- bg_direct_acceptance(key_summary, motif_detail)
    benchmark <- NULL
    if (!is.null(benchmark_dir) && nzchar(benchmark_dir)) {
        benchmark_key <- pairs$reference$organism == "hs" &
            pairs$reference$assembly == "hg38" &
            pairs$reference$upstream == 950L &
            pairs$reference$downstream == 50L &
            pairs$reference$jaspar_release == 2020L
        if (sum(benchmark_key) != 1L) {
            bg_stop("Cannot resolve the hg38 JASPAR2020 benchmark background")
        }
        reference_bg <- bg_read_background_table(file.path(
            root, pairs$reference$artifact[benchmark_key]
        ))
        candidate_bg <- bg_read_background_table(file.path(
            root, pairs$candidate$artifact[benchmark_key]
        ))
        benchmark <- bg_compare_benchmark(
            benchmark_dir, reference_bg, candidate_bg
        )
        checks <- rbind(checks, bg_benchmark_acceptance(benchmark))
    }
    result <- list(
        key_summary = key_summary, motif_detail = motif_detail,
        benchmark = benchmark, checks = checks
    )
    bg_write_comparison_reports(
        result, reference_version, candidate_version, output_dir
    )
    result
}
