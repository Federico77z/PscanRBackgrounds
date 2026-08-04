BG_SCHEMA_VERSION <- 1L
BG_SCORING_SPEC <- "pscan-best-hit-v1"
BG_CATALOG_COLUMNS <- c(
    "schema_version", "status", "latest", "organism", "assembly",
    "upstream", "downstream", "jaspar_release", "tax_group",
    "background_version", "artifact", "artifact_sha256", "motif_count",
    "promoter_count", "annotation_hash", "promoter_hash", "motif_hash",
    "annotation_snapshot", "annotation_source", "annotation_retrieved_at",
    "scoring_spec", "pscanr_version", "pscanr_commit", "bioc_version",
    "generated_at"
)

bg_stop <- function(...) stop(..., call. = FALSE)

bg_message <- function(...) {
    message(sprintf(...))
}

bg_require <- function(packages) {
    missing <- packages[!vapply(
        packages, requireNamespace, logical(1), quietly = TRUE
    )]
    if (length(missing)) {
        bg_stop("Missing required R packages: ", paste(missing, collapse = ", "))
    }
}

bg_read_config <- function(root) {
    read_tsv <- function(name) {
        utils::read.delim(
            file.path(root, "config", name), stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
    list(
        organisms = read_tsv("backgrounds.tsv"),
        windows = read_tsv("windows.tsv"),
        jaspar = read_tsv("jaspar.tsv")
    )
}

bg_empty_catalog <- function() {
    out <- as.data.frame(
        setNames(replicate(length(BG_CATALOG_COLUMNS), character(), FALSE),
                 BG_CATALOG_COLUMNS),
        stringsAsFactors = FALSE
    )
    for (column in c(
        "schema_version", "upstream", "downstream", "jaspar_release",
        "background_version", "motif_count", "promoter_count"
    )) {
        out[[column]] <- integer()
    }
    out$latest <- logical()
    out
}

bg_read_catalog <- function(root) {
    path <- file.path(root, "catalog.tsv")
    if (!file.exists(path)) return(bg_empty_catalog())
    out <- utils::read.delim(
        path, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character()
    )
    missing <- setdiff(BG_CATALOG_COLUMNS, names(out))
    if (length(missing)) {
        bg_stop("catalog.tsv is missing columns: ", paste(missing, collapse = ", "))
    }
    out <- out[, BG_CATALOG_COLUMNS, drop = FALSE]
    integer_columns <- c(
        "schema_version", "upstream", "downstream", "jaspar_release",
        "background_version", "motif_count", "promoter_count"
    )
    out[integer_columns] <- lapply(out[integer_columns], as.integer)
    out$latest <- tolower(out$latest) == "true"
    out
}

bg_write_catalog <- function(catalog, root) {
    catalog <- catalog[order(
        catalog$jaspar_release, catalog$assembly, catalog$upstream,
        catalog$downstream, catalog$background_version
    ), BG_CATALOG_COLUMNS, drop = FALSE]
    path <- file.path(root, "catalog.tsv")
    temporary <- paste0(path, ".tmp")
    utils::write.table(
        catalog, temporary, sep = "\t", quote = FALSE,
        row.names = FALSE, na = ""
    )
    if (!file.rename(temporary, path)) bg_stop("Could not replace ", path)

    available <- catalog$artifact[
        catalog$status == "validated" & catalog$latest
    ]
    writeLines(sort(basename(available)), file.path(root, "AvailableBG.txt"))
    invisible(path)
}

bg_hash_file <- function(path) {
    if (!file.exists(path)) bg_stop("Cannot hash missing file: ", path)
    unname(tools::sha256sum(path))
}

bg_hash_object <- function(object) {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(object, path, version = 3, compress = "xz")
    bg_hash_file(path)
}

bg_key <- function(organism, assembly, upstream, downstream, jaspar_release) {
    paste(organism, assembly, upstream, downstream, jaspar_release, sep = ":")
}

bg_filename <- function(jaspar_release, assembly, upstream, downstream,
                        provider, version) {
    suffix <- if (identical(provider, "gff") && identical(assembly, "TAIR9")) {
        "TAIR"
    } else {
        "UCSC"
    }
    sprintf(
        "J%d_%s_%du_%dd_%s.psbg%d.txt", jaspar_release, assembly,
        upstream, downstream, suffix, version
    )
}

bg_filter_rows <- function(config, filters) {
    organisms <- config$organisms
    windows <- config$windows
    jaspar <- config$jaspar
    if (length(filters$assembly)) {
        organisms <- organisms[organisms$assembly %in% filters$assembly, ]
    }
    if (length(filters$jaspar)) {
        jaspar <- jaspar[jaspar$release %in% as.integer(filters$jaspar), ]
    }
    if (length(filters$window)) {
        labels <- paste0(windows$upstream, "u_", windows$downstream, "d")
        windows <- windows[labels %in% filters$window, ]
    }
    if (!nrow(organisms) || !nrow(windows) || !nrow(jaspar)) {
        bg_stop("Filters selected no configured background combinations")
    }
    list(organisms = organisms, windows = windows, jaspar = jaspar)
}

bg_package_object <- function(package) {
    bg_require(package)
    getExportedValue(package, package)
}

bg_name_transcripts <- function(tx, curated_refseq = FALSE) {
    metadata <- S4Vectors::mcols(tx)
    tx_name <- as.character(metadata$tx_name)
    missing_name <- is.na(tx_name) | !nzchar(tx_name)
    if (curated_refseq) {
        keep <- !missing_name & !grepl("^X[MR]_", tx_name)
        tx <- tx[keep]
        names(tx) <- tx_name[keep]
        return(tx)
    }
    tx_id <- as.character(metadata$tx_id)
    if (length(tx_id) != length(tx)) tx_id <- rep(NA_character_, length(tx))
    usable_id <- !is.na(tx_id) & nzchar(tx_id)
    tx_name[missing_name & usable_id] <- paste0("tx-", tx_id[missing_name & usable_id])
    missing_name <- is.na(tx_name) | !nzchar(tx_name)
    tx_name[missing_name] <- paste0("tx-", which(missing_name))
    names(tx) <- tx_name
    tx
}

bg_load_annotation <- function(spec, root) {
    bg_require(c("GenomicFeatures", "GenomeInfoDb", "txdbmaker"))
    provider <- spec$provider[[1]]
    source <- spec$annotation_source[[1]]
    assembly <- spec$assembly[[1]]

    if (provider == "ucsc") {
        txdb <- txdbmaker::makeTxDbFromUCSC(
            genome = assembly, tablename = source
        )
    } else if (provider %in% c("gff", "ncbi_gff")) {
        bg_require("rtracklayer")
        path <- file.path(root, source)
        if (!file.exists(path)) bg_stop("Missing annotation source: ", path)
        if (provider == "ncbi_gff") {
            gff <- rtracklayer::import(path)
            regions <- gff[as.character(S4Vectors::mcols(gff)$type) == "region"]
            chromosome <- as.character(S4Vectors::mcols(regions)$chromosome)
            mapping <- stats::setNames(
                paste0("chr", chromosome),
                as.character(GenomicRanges::seqnames(regions))
            )
            mapping <- mapping[!is.na(chromosome) & nzchar(chromosome)]
            present <- intersect(names(mapping), GenomeInfoDb::seqlevels(gff))
            GenomeInfoDb::seqlevels(gff) <- mapping[present]
            txdb <- txdbmaker::makeTxDbFromGRanges(gff, taxonomyId = 9606)
        } else {
            txdb <- txdbmaker::makeTxDbFromGFF(
                path, format = "gff3", dataSource = assembly,
                organism = "Arabidopsis thaliana"
            )
        }
    } else {
        bg_stop("Unsupported annotation provider: ", provider)
    }

    tx <- GenomicFeatures::transcripts(
        txdb, columns = c("tx_id", "tx_name", "gene_id")
    )
    canonical <- strsplit(
        spec$canonical_chromosomes[[1]], ",", fixed = TRUE
    )[[1]]
    tx <- tx[as.character(GenomicRanges::seqnames(tx)) %in% canonical]
    tx <- bg_name_transcripts(
        tx, spec$transcript_filter[[1]] == "curated_refseq"
    )
    if (!length(tx)) bg_stop("No transcripts retained for ", assembly)
    tx
}

bg_annotation_table <- function(tx) {
    genes <- S4Vectors::mcols(tx)$gene_id
    gene_id <- vapply(genes, function(x) {
        if (!length(x)) "" else paste(sort(as.character(x)), collapse = ",")
    }, character(1))
    strand <- as.character(BiocGenerics::strand(tx))
    tss <- ifelse(strand == "+", BiocGenerics::start(tx), BiocGenerics::end(tx))
    data.frame(
        transcript_id = names(tx), gene_id = gene_id,
        seqname = as.character(GenomicRanges::seqnames(tx)),
        strand = strand, tss = as.integer(tss),
        stringsAsFactors = FALSE
    )
}

bg_annotation_hash <- function(tx) {
    tab <- unique(bg_annotation_table(tx)[, c("seqname", "strand", "tss")])
    tab <- tab[order(tab$seqname, tab$strand, tab$tss), ]
    row.names(tab) <- NULL
    bg_hash_object(tab)
}

bg_prepare_promoters <- function(spec, upstream, downstream, root,
                                 write_cache = TRUE, tx = NULL) {
    bg_require(c("Biostrings", "BSgenome", "GenomicRanges"))
    if (is.null(tx)) tx <- bg_load_annotation(spec, root)
    annotation_hash <- bg_annotation_hash(tx)
    genome_package <- spec$bsgenome[[1]]
    genome <- bg_package_object(genome_package)
    ranges <- GenomicRanges::promoters(
        tx, upstream = upstream, downstream = downstream
    )
    lengths <- GenomeInfoDb::seqlengths(genome)
    seqnames <- as.character(GenomicRanges::seqnames(ranges))
    chromosome_lengths <- unname(lengths[seqnames])
    valid <- !is.na(chromosome_lengths) & BiocGenerics::start(ranges) >= 1L &
        BiocGenerics::end(ranges) <= chromosome_lengths
    excluded_out_of_bounds <- sum(!valid)
    ranges <- ranges[valid]
    tx <- tx[valid]
    sequences <- Biostrings::getSeq(genome, ranges)
    names(sequences) <- names(tx)
    expected_width <- upstream + downstream
    exact_width <- Biostrings::width(sequences) == expected_width
    excluded_wrong_width <- sum(!exact_width)
    sequences <- sequences[exact_width]
    tx <- tx[exact_width]
    n_fraction <- Biostrings::letterFrequency(sequences, "N", as.prob = TRUE)
    acceptable <- n_fraction <= 0.5
    excluded_high_n <- sum(!acceptable)
    sequences <- sequences[acceptable]
    tx <- tx[acceptable]
    raw_count <- length(sequences)
    all_sequences <- sequences
    sequences <- BiocGenerics::unique(sequences)
    unique_count <- length(sequences)
    sorted_sequences <- sort(as.character(sequences))
    promoter_hash <- bg_hash_object(sorted_sequences)

    annotation <- bg_annotation_table(tx)
    snapshot_name <- sprintf(
        "%s_%du_%dd_%s.rds", spec$assembly[[1]], upstream, downstream,
        substr(promoter_hash, 1L, 16L)
    )
    snapshot <- list(
        schema_version = BG_SCHEMA_VERSION,
        assembly = spec$assembly[[1]], source = spec$annotation_source[[1]],
        annotation_hash = annotation_hash, promoter_hash = promoter_hash,
        upstream = upstream, downstream = downstream,
        transcript_count = raw_count, unique_promoter_count = unique_count,
        exclusions = c(out_of_bounds = excluded_out_of_bounds,
                       wrong_width = excluded_wrong_width,
                       high_n = excluded_high_n),
        transcripts = annotation
    )
    cache_path <- file.path(
        root, ".cache", paste0("promoters_", promoter_hash, ".rds")
    )
    if (write_cache && !file.exists(cache_path)) {
        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(
            list(
                sequences = sequences, all_sequences = all_sequences,
                snapshot = snapshot
            ),
            cache_path, compress = "xz"
        )
    }
    list(
        sequences = sequences, snapshot = snapshot,
        snapshot_name = snapshot_name, cache_path = cache_path,
        annotation_hash = annotation_hash, promoter_hash = promoter_hash,
        promoter_count = unique_count
    )
}

bg_load_motifs <- function(release, tax_group) {
    package <- paste0("JASPAR", release)
    bg_require(c(package, "TFBSTools"))
    opts <- list(collection = "CORE", tax_group = tax_group)
    if (release == 2024L) {
        bg_require("RSQLite")
        database <- getExportedValue(package, package)()
        connection <- RSQLite::dbConnect(
            RSQLite::SQLite(), getExportedValue(package, "db")(database)
        )
        on.exit(RSQLite::dbDisconnect(connection), add = TRUE)
        motifs <- TFBSTools::getMatrixSet(connection, opts)
    } else {
        motifs <- TFBSTools::getMatrixSet(bg_package_object(package), opts)
    }
    ids <- TFBSTools::ID(motifs)
    matrices <- lapply(motifs, function(x) unclass(TFBSTools::Matrix(x)))
    names(matrices) <- ids
    list(
        motifs = motifs, motif_hash = bg_hash_object(matrices),
        motif_count = length(motifs)
    )
}

bg_latest_for_key <- function(catalog, key) {
    if (!nrow(catalog)) return(catalog)
    keys <- bg_key(
        catalog$organism, catalog$assembly, catalog$upstream,
        catalog$downstream, catalog$jaspar_release
    )
    catalog[keys == key & catalog$latest & catalog$status == "validated", ]
}

bg_next_version <- function(catalog, key) {
    if (!nrow(catalog)) return(1L)
    keys <- bg_key(
        catalog$organism, catalog$assembly, catalog$upstream,
        catalog$downstream, catalog$jaspar_release
    )
    versions <- catalog$background_version[keys == key]
    if (!length(versions)) 1L else max(versions, na.rm = TRUE) + 1L
}

bg_build_plan <- function(config, catalog, root, filters, force = FALSE) {
    selected <- bg_filter_rows(config, filters)
    rows <- list()
    index <- 0L
    motif_cache <- new.env(parent = emptyenv())
    for (i in seq_len(nrow(selected$organisms))) {
        spec <- selected$organisms[i, , drop = FALSE]
        tx <- bg_load_annotation(spec, root)
        for (j in seq_len(nrow(selected$windows))) {
            upstream <- selected$windows$upstream[[j]]
            downstream <- selected$windows$downstream[[j]]
            prepared <- bg_prepare_promoters(
                spec, upstream, downstream, root, write_cache = TRUE, tx = tx
            )
            for (k in seq_len(nrow(selected$jaspar))) {
                release <- selected$jaspar$release[[k]]
                motif_key <- paste(release, spec$tax_group[[1]], sep = ":")
                if (!exists(motif_key, motif_cache, inherits = FALSE)) {
                    assign(
                        motif_key,
                        bg_load_motifs(release, spec$tax_group[[1]]),
                        motif_cache
                    )
                }
                motif <- get(motif_key, motif_cache, inherits = FALSE)
                key <- bg_key(
                    spec$organism[[1]], spec$assembly[[1]], upstream,
                    downstream, release
                )
                current <- bg_latest_for_key(catalog, key)
                unchanged <- nrow(current) == 1L &&
                    nzchar(current$promoter_hash) &&
                    identical(current$promoter_hash, prepared$promoter_hash) &&
                    identical(current$motif_hash, motif$motif_hash) &&
                    identical(current$scoring_spec, BG_SCORING_SPEC)
                index <- index + 1L
                rows[[index]] <- data.frame(
                    key = key, action = if (unchanged && !force) "skip" else "generate",
                    reason = if (force) "forced" else if (unchanged) "unchanged" else
                        if (!nrow(current) || !nzchar(current$promoter_hash))
                            "legacy_or_missing_fingerprint" else "input_changed",
                    organism = spec$organism[[1]], assembly = spec$assembly[[1]],
                    provider = spec$provider[[1]], tax_group = spec$tax_group[[1]],
                    annotation_source = spec$annotation_source[[1]],
                    upstream = upstream, downstream = downstream,
                    jaspar_release = release,
                    background_version = bg_next_version(catalog, key),
                    promoter_hash = prepared$promoter_hash,
                    annotation_hash = prepared$annotation_hash,
                    promoter_count = prepared$promoter_count,
                    promoter_cache = prepared$cache_path,
                    snapshot_name = prepared$snapshot_name,
                    motif_hash = motif$motif_hash,
                    motif_count = motif$motif_count,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    do.call(rbind, rows)
}

bg_validate_artifact <- function(path, expected_motif_count = NA_integer_) {
    problems <- character()
    if (!file.exists(path)) return("file is missing")
    first <- readLines(path, n = 1L, warn = FALSE)
    if (!identical(first, "[SHORT TFBS MATRIX]")) {
        problems <- c(problems, "invalid header")
    }
    table <- tryCatch(
        utils::read.delim(path, skip = 1L, header = FALSE),
        error = function(e) e
    )
    if (inherits(table, "error")) return(conditionMessage(table))
    if (ncol(table) != 4L) problems <- c(problems, "expected four columns")
    if (!is.na(expected_motif_count) && nrow(table) != expected_motif_count) {
        problems <- c(problems, sprintf(
            "expected %d motifs, found %d", expected_motif_count, nrow(table)
        ))
    }
    if (ncol(table) == 4L && nrow(table)) {
        if (anyDuplicated(table[[1]])) problems <- c(problems, "duplicate motif IDs")
        if (length(unique(table[[2]])) != 1L || any(table[[2]] <= 0)) {
            problems <- c(problems, "invalid background size")
        }
        if (any(!is.finite(table[[3]]) | table[[3]] < 0 | table[[3]] > 1)) {
            problems <- c(problems, "invalid background mean")
        }
        if (any(!is.finite(table[[4]]) | table[[4]] <= 0 | table[[4]] > 1)) {
            problems <- c(problems, "invalid background standard deviation")
        }
    }
    unique(problems)
}

bg_git_value <- function(root, args) {
    value <- tryCatch(
        system2("git", c("-C", shQuote(root), args), stdout = TRUE, stderr = FALSE),
        error = function(e) character()
    )
    if (length(value)) value[[1]] else ""
}

bg_pscanr_metadata <- function() {
    source <- Sys.getenv("PSCANR_SOURCE", "")
    if (nzchar(source) && file.exists(file.path(source, "DESCRIPTION"))) {
        description <- read.dcf(file.path(source, "DESCRIPTION"))
        list(
            version = unname(description[1, "Version"]),
            commit = bg_git_value(source, c("rev-parse", "HEAD"))
        )
    } else if (requireNamespace("PscanR", quietly = TRUE)) {
        list(version = as.character(utils::packageVersion("PscanR")), commit = "")
    } else {
        bg_stop("Install PscanR or set PSCANR_SOURCE to a clean PscanR checkout")
    }
}

bg_load_pscanr <- function() {
    source <- Sys.getenv("PSCANR_SOURCE", "")
    if (nzchar(source)) {
        bg_require("pkgload")
        status <- system2(
            "git", c("-C", shQuote(source), "status", "--porcelain"),
            stdout = TRUE
        )
        if (length(status)) bg_stop("PSCANR_SOURCE must be a clean checkout")
        pkgload::load_all(source, quiet = TRUE, export_all = FALSE)
    } else {
        bg_require("PscanR")
    }
    invisible(TRUE)
}

bg_publish_snapshot <- function(plan_row, prepared, root) {
    relative <- file.path("annotation_snapshots", plan_row$snapshot_name)
    path <- file.path(root, relative)
    if (!file.exists(path)) {
        saveRDS(prepared$snapshot, path, version = 3, compress = "xz")
    }
    relative
}

bg_generate <- function(plan, catalog, root, cores) {
    jobs <- plan[plan$action == "generate", , drop = FALSE]
    if (!nrow(jobs)) {
        bg_message("No background regeneration is required")
        return(catalog)
    }
    bg_load_pscanr()
    metadata <- bg_pscanr_metadata()
    bpparam <- if (cores == 1L) {
        BiocParallel::SerialParam()
    } else if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(cores)
    } else {
        BiocParallel::MulticoreParam(cores)
    }

    for (i in seq_len(nrow(jobs))) {
        job <- jobs[i, , drop = FALSE]
        started_at <- Sys.time()
        spec <- bg_read_config(root)$organisms
        spec <- spec[spec$assembly == job$assembly, , drop = FALSE]
        prepared <- readRDS(job$promoter_cache)
        motif <- bg_load_motifs(job$jaspar_release, job$tax_group)
        filename <- bg_filename(
            job$jaspar_release, job$assembly, job$upstream, job$downstream,
            spec$provider[[1]], job$background_version
        )
        stage <- file.path(root, "staging", filename)
        final <- file.path(root, "BG_files", filename)
        bg_message("[%d/%d] Generating %s", i, nrow(jobs), filename)
        result <- PscanR::ps_build_bg(
            prepared$sequences, motif$motifs, BPPARAM = bpparam
        )
        PscanR::ps_write_bg_to_file(result, stage)
        problems <- bg_validate_artifact(stage, job$motif_count)
        if (length(problems)) {
            bg_stop("Validation failed for ", filename, ": ",
                    paste(problems, collapse = "; "))
        }
        if (file.exists(final)) bg_stop("Refusing to overwrite ", final)
        if (!file.rename(stage, final)) bg_stop("Could not publish ", filename)
        snapshot <- bg_publish_snapshot(job, prepared, root)

        key <- job$key
        if (nrow(catalog)) {
            keys <- bg_key(
                catalog$organism, catalog$assembly, catalog$upstream,
                catalog$downstream, catalog$jaspar_release
            )
            catalog$latest[keys == key] <- FALSE
        }
        new <- as.list(setNames(rep("", length(BG_CATALOG_COLUMNS)),
                                BG_CATALOG_COLUMNS))
        new$schema_version <- BG_SCHEMA_VERSION
        new$status <- "validated"
        new$latest <- TRUE
        for (field in c(
            "organism", "assembly", "upstream", "downstream",
            "jaspar_release", "tax_group", "background_version",
            "motif_count", "promoter_count", "annotation_hash",
            "promoter_hash", "motif_hash"
        )) new[[field]] <- job[[field]][[1]]
        new$artifact <- file.path("BG_files", filename)
        new$artifact_sha256 <- bg_hash_file(final)
        new$annotation_snapshot <- snapshot
        new$annotation_source <- job$annotation_source
        new$annotation_retrieved_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
        new$scoring_spec <- BG_SCORING_SPEC
        new$pscanr_version <- metadata$version
        new$pscanr_commit <- metadata$commit
        new$bioc_version <- as.character(BiocManager::version())
        new$generated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
        catalog <- rbind(catalog, as.data.frame(new, stringsAsFactors = FALSE))
        bg_write_catalog(catalog, root)
        bg_message(
            "[%d/%d] Completed %s in %.1f seconds",
            i, nrow(jobs), filename,
            as.numeric(difftime(Sys.time(), started_at, units = "secs"))
        )
        rm(result)
        gc()
    }
    catalog
}

bg_validate_catalog <- function(catalog, root) {
    problems <- list()
    seen_keys <- character()
    for (i in seq_len(nrow(catalog))) {
        row <- catalog[i, ]
        path <- file.path(root, row$artifact)
        issue <- bg_validate_artifact(path, row$motif_count)
        if (file.exists(path) && nzchar(row$artifact_sha256) &&
            !identical(bg_hash_file(path), row$artifact_sha256)) {
            issue <- c(issue, "SHA-256 mismatch")
        }
        key <- paste(
            bg_key(row$organism, row$assembly, row$upstream, row$downstream,
                   row$jaspar_release), row$background_version, sep = ":"
        )
        if (key %in% seen_keys) issue <- c(issue, "duplicate catalog key/version")
        seen_keys <- c(seen_keys, key)
        if (length(issue)) problems[[row$artifact]] <- unique(issue)
    }
    latest <- catalog[catalog$latest & catalog$status == "validated", ]
    latest_keys <- bg_key(
        latest$organism, latest$assembly, latest$upstream,
        latest$downstream, latest$jaspar_release
    )
    duplicate_latest <- unique(latest_keys[duplicated(latest_keys)])
    if (length(duplicate_latest)) {
        problems[["catalog.tsv"]] <- paste(
            "multiple latest entries:", paste(duplicate_latest, collapse = ", ")
        )
    }
    problems
}

bg_audit_existing <- function(root, catalog) {
    files <- list.files(file.path(root, "BG_files"), full.names = TRUE)
    rows <- lapply(files, function(path) {
        table <- utils::read.delim(path, skip = 1L, header = FALSE)
        data.frame(
            artifact = basename(path), motif_count = nrow(table),
            promoter_count = if (nrow(table)) table[[2]][[1]] else NA_integer_,
            problems = paste(bg_validate_artifact(path), collapse = ";"),
            sha256 = bg_hash_file(path), stringsAsFactors = FALSE
        )
    })
    audit <- do.call(rbind, rows)
    utils::write.table(
        audit, file.path(root, "reports", "existing_artifact_audit.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    catalog_problems <- bg_validate_catalog(catalog, root)
    bg_message(
        "Audited %d files: %d artifact failures, %d catalog failures",
        nrow(audit), sum(nzchar(audit$problems)), length(catalog_problems)
    )
    invisible(list(artifacts = audit, catalog = catalog_problems))
}

bg_check_annotations <- function(config, catalog, root, filters) {
    selected <- bg_filter_rows(config, filters)
    rows <- lapply(seq_len(nrow(selected$organisms)), function(i) {
        spec <- selected$organisms[i, , drop = FALSE]
        tx <- bg_load_annotation(spec, root)
        hash <- bg_annotation_hash(tx)
        previous <- unique(catalog$annotation_hash[
            catalog$assembly == spec$assembly & catalog$latest
        ])
        previous <- previous[nzchar(previous)]
        data.frame(
            assembly = spec$assembly, annotation_source = spec$annotation_source,
            current_annotation_hash = hash,
            previous_annotation_hash = paste(previous, collapse = ","),
            candidate_change = !length(previous) || !hash %in% previous,
            checked_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    utils::write.table(
        out, file.path(root, "reports", "annotation_check.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    print(out, row.names = FALSE)
    invisible(out)
}

bg_stratified_motif_indices <- function(motifs, count = 50L) {
    widths <- vapply(motifs, function(x) ncol(TFBSTools::Matrix(x)), integer(1))
    information <- vapply(motifs, function(x) {
        matrix <- TFBSTools::Matrix(x)
        probability <- sweep(matrix, 2L, colSums(matrix), "/")
        terms <- ifelse(probability > 0, probability * log2(probability / 0.25), 0)
        sum(terms)
    }, numeric(1))
    bins <- interaction(
        cut(rank(widths, ties.method = "first"), 5L, labels = FALSE),
        cut(rank(information, ties.method = "first"), 5L, labels = FALSE),
        drop = TRUE
    )
    selected <- unlist(lapply(split(seq_along(motifs), bins), head, 2L))
    if (length(selected) < min(count, length(motifs))) {
        remaining <- setdiff(seq_along(motifs), selected)
        needed <- min(count, length(motifs)) - length(selected)
        positions <- unique(round(seq(1, length(remaining), length.out = needed)))
        selected <- c(selected, remaining[positions])
    }
    sort(head(unique(selected), count))
}

bg_scan_score_matrix <- function(sequences, motifs, cores) {
    parameter <- if (cores == 1L) {
        BiocParallel::SerialParam()
    } else if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(cores)
    } else {
        BiocParallel::MulticoreParam(cores)
    }
    scores <- BiocParallel::bplapply(
        motifs,
        function(motif, sequences) {
            psm <- methods::as(motif, "PSMatrix")
            hit <- PscanR::ps_scan(
                psm, sequences, BG = TRUE, fullBG = TRUE
            )
            PscanR::ps_hits_score_bg(hit)
        },
        sequences = sequences, BPPARAM = parameter
    )
    matrix(
        unlist(scores, use.names = FALSE), nrow = length(sequences),
        ncol = length(motifs), dimnames = list(names(sequences), names(motifs))
    )
}

bg_one_promoter_per_gene <- function(prepared) {
    sequences <- prepared$all_sequences
    annotation <- prepared$snapshot$transcripts
    annotation <- annotation[match(names(sequences), annotation$transcript_id), ]
    gene <- sub(",.*$", "", annotation$gene_id)
    gene[!nzchar(gene)] <- annotation$transcript_id[!nzchar(gene)]
    priority <- ifelse(
        grepl("^NM_", annotation$transcript_id), 1L,
        ifelse(grepl("^NR_", annotation$transcript_id), 2L, 3L)
    )
    ordering <- order(gene, priority, annotation$transcript_id)
    keep <- ordering[!duplicated(gene[ordering])]
    BiocGenerics::unique(sequences[keep])
}

bg_calibrate <- function(config, root, filters, cores, repetitions = 1000L,
                         seed = 20090516L) {
    bg_load_pscanr()
    selected <- bg_filter_rows(config, filters)$organisms
    details <- list()
    comparisons <- list()
    detail_index <- 0L
    comparison_index <- 0L
    set_sizes <- c(5L, 10L, 20L, 50L, 100L, 200L)

    for (i in seq_len(nrow(selected))) {
        spec <- selected[i, , drop = FALSE]
        assembly <- spec$assembly[[1]]
        bg_message("Scientific calibration for %s", assembly)
        prepared_info <- bg_prepare_promoters(
            spec, 450L, 50L, root, write_cache = TRUE
        )
        prepared <- readRDS(prepared_info$cache_path)
        motif_info <- bg_load_motifs(2024L, spec$tax_group[[1]])
        scores <- bg_scan_score_matrix(
            prepared$sequences, motif_info$motifs, cores
        )
        means <- colMeans(scores, na.rm = TRUE)
        deviations <- apply(scores, 2L, stats::sd, na.rm = TRUE)
        set.seed(seed + i)

        for (set_size in set_sizes[set_sizes <= nrow(scores)]) {
            samples <- replicate(
                repetitions, sample.int(nrow(scores), set_size),
                simplify = FALSE
            )
            for (motif_index in seq_len(ncol(scores))) {
                values <- scores[, motif_index]
                sample_means <- vapply(
                    samples, function(index) mean(values[index], na.rm = TRUE),
                    numeric(1)
                )
                z <- (sample_means - means[[motif_index]]) /
                    (deviations[[motif_index]] / sqrt(set_size))
                pvalues <- stats::pnorm(z, lower.tail = FALSE)
                pvalues <- pvalues[is.finite(pvalues)]
                detail_index <- detail_index + 1L
                details[[detail_index]] <- data.frame(
                    assembly = assembly,
                    motif_id = names(motif_info$motifs)[[motif_index]],
                    set_size = set_size,
                    ks_pvalue = if (length(pvalues) >= 10L) {
                        suppressWarnings(stats::ks.test(pvalues, "punif")$p.value)
                    } else {
                        NA_real_
                    },
                    fpr_005 = mean(pvalues < 0.05),
                    fpr_001 = mean(pvalues < 0.01),
                    repetitions = repetitions,
                    stringsAsFactors = FALSE
                )
            }
        }

        diagnostic_indices <- bg_stratified_motif_indices(motif_info$motifs)
        diagnostic_motifs <- motif_info$motifs[diagnostic_indices]
        gene_sequences <- bg_one_promoter_per_gene(prepared)
        gene_scores <- bg_scan_score_matrix(gene_sequences, diagnostic_motifs, cores)
        for (j in seq_along(diagnostic_indices)) {
            current_index <- diagnostic_indices[[j]]
            comparison_index <- comparison_index + 1L
            comparisons[[comparison_index]] <- data.frame(
                assembly = assembly,
                motif_id = names(motif_info$motifs)[[current_index]],
                transcript_promoters = nrow(scores),
                gene_promoters = nrow(gene_scores),
                transcript_mean = means[[current_index]],
                gene_mean = mean(gene_scores[, j], na.rm = TRUE),
                transcript_sd = deviations[[current_index]],
                gene_sd = stats::sd(gene_scores[, j], na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        }
        rm(scores, gene_scores)
        gc()
    }

    details <- do.call(rbind, details)
    details$ks_fdr <- stats::p.adjust(details$ks_pvalue, method = "BH")
    comparisons <- do.call(rbind, comparisons)
    comparisons$mean_delta <- comparisons$gene_mean - comparisons$transcript_mean
    comparisons$sd_delta <- comparisons$gene_sd - comparisons$transcript_sd
    utils::write.table(
        details, file.path(root, "reports", "scientific_calibration.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
        comparisons, file.path(root, "reports", "promoter_universe_comparison.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    grouped <- split(details, interaction(details$assembly, details$set_size))
    summary <- do.call(rbind, lapply(grouped, function(x) {
        n <- nrow(x) * x$repetitions[[1]]
        limit_005 <- 0.05 + 3 * sqrt(0.05 * 0.95 / n)
        limit_001 <- 0.01 + 3 * sqrt(0.01 * 0.99 / n)
        data.frame(
            assembly = x$assembly[[1]], set_size = x$set_size[[1]],
            motifs = nrow(x),
            rejected_fraction = mean(x$ks_fdr < 0.05, na.rm = TRUE),
            aggregate_fpr_005 = mean(x$fpr_005, na.rm = TRUE),
            aggregate_fpr_001 = mean(x$fpr_001, na.rm = TRUE),
            limit_005 = limit_005, limit_001 = limit_001,
            passed = mean(x$ks_fdr < 0.05, na.rm = TRUE) <= 0.05 &&
                mean(x$fpr_005, na.rm = TRUE) <= limit_005 &&
                mean(x$fpr_001, na.rm = TRUE) <= limit_001,
            stringsAsFactors = FALSE
        )
    }))
    row.names(summary) <- NULL
    utils::write.table(
        summary, file.path(root, "reports", "scientific_calibration_summary.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    list(details = details, comparison = comparisons, summary = summary)
}
