BG_ZENODO_ARCHIVE_DATE <- as.POSIXct(
    "2020-01-01 00:00:00", tz = "UTC"
)

bg_zenodo_release_catalog <- function(catalog, version) {
    selected <- catalog[
        catalog$status == "validated" &
            catalog$background_version == version,
        , drop = FALSE
    ]
    if (!nrow(selected)) {
        bg_stop("No validated backgrounds found for version ", version)
    }
    keys <- bg_key(
        selected$organism, selected$assembly, selected$upstream,
        selected$downstream, selected$jaspar_release
    )
    if (anyDuplicated(keys)) {
        bg_stop("Background version ", version, " has duplicate catalog keys")
    }
    selected <- selected[order(
        selected$jaspar_release, selected$assembly, selected$upstream,
        selected$downstream
    ), BG_CATALOG_COLUMNS, drop = FALSE]
    selected$artifact <- file.path("backgrounds", basename(selected$artifact))
    row.names(selected) <- NULL
    selected
}

bg_zenodo_readme <- function(version, catalog) {
    paste(
        "# PscanR precomputed promoter backgrounds",
        "",
        sprintf("This archive contains immutable PscanR background version %d.", version),
        sprintf(
            "It provides %d validated combinations of JASPAR release, genome assembly, and promoter window.",
            nrow(catalog)
        ),
        "",
        "The files store the mean and standard deviation of the best normalized PWM score across the unique promoter sequences used for each background.",
        "They are intended for use by the PscanR motif-enrichment package.",
        "",
        "## Contents",
        "",
        "- `backgrounds/`: PscanR short-background text files.",
        "- `catalog.tsv`: resource keys, provenance, and SHA-256 checksums.",
        "- `MANIFEST.sha256`: SHA-256 checksums for all other archive members.",
        "- `CITATION.cff`: citation metadata.",
        "- `LICENSE`: CC BY 4.0 license notice.",
        "",
        "## Supported resources",
        "",
        paste0("- JASPAR releases: ", paste(sort(unique(catalog$jaspar_release)), collapse = ", "), "."),
        paste0("- Assemblies: ", paste(sort(unique(catalog$assembly)), collapse = ", "), "."),
        paste0(
            "- Promoter windows: ",
            paste(unique(sprintf("%du_%dd", catalog$upstream, catalog$downstream)), collapse = ", "),
            "."
        ),
        "",
        "## Format",
        "",
        "Each background starts with `[SHORT TFBS MATRIX]`, followed by tab-separated motif identifier, promoter count, mean score, and score standard deviation.",
        "The `artifact_sha256` column in `catalog.tsv` validates each background file.",
        "",
        "Generation and scientific-validation code is maintained at:",
        "https://github.com/Federico77z/PscanR_backgrounds",
        "",
        "PscanR is maintained at:",
        "https://github.com/Federico77z/PscanR",
        "",
        "Please cite Zambelli F, Pesole G, Pavesi G. Pscan: finding over-represented transcription factor binding site motifs in sequences from co-regulated or co-expressed genes. Nucleic Acids Research (2009). https://doi.org/10.1093/nar/gkp464",
        sep = "\n"
    )
}

bg_zenodo_license <- function() {
    paste(
        "PscanR precomputed promoter backgrounds",
        "",
        "Copyright (c) Federico Zambelli and Giulio Pavesi.",
        "",
        "This dataset is licensed under the Creative Commons Attribution 4.0 International license (CC BY 4.0).",
        "",
        "You are free to share and adapt the material for any purpose, provided appropriate credit is given, a link to the license is supplied, and changes are indicated.",
        "",
        "License text: https://creativecommons.org/licenses/by/4.0/legalcode",
        sep = "\n"
    )
}

bg_zenodo_citation <- function(version) {
    paste(
        "cff-version: 1.2.0",
        "message: \"If you use these backgrounds, please cite the Pscan publication and this dataset.\"",
        "title: \"PscanR precomputed promoter backgrounds\"",
        sprintf("version: \"%d\"", version),
        "type: dataset",
        "authors:",
        "  - family-names: Zambelli",
        "    given-names: Federico",
        "    orcid: \"https://orcid.org/0000-0003-3487-4331\"",
        "  - family-names: Pavesi",
        "    given-names: Giulio",
        "    orcid: \"https://orcid.org/0000-0001-5705-6249\"",
        "license: CC-BY-4.0",
        "repository-code: \"https://github.com/Federico77z/PscanR_backgrounds\"",
        "references:",
        "  - type: article",
        "    authors:",
        "      - family-names: Zambelli",
        "        given-names: Federico",
        "      - family-names: Pesole",
        "        given-names: Graziano",
        "      - family-names: Pavesi",
        "        given-names: Giulio",
        "    title: \"Pscan: finding over-represented transcription factor binding site motifs in sequences from co-regulated or co-expressed genes\"",
        "    journal: \"Nucleic Acids Research\"",
        "    year: 2009",
        "    doi: \"10.1093/nar/gkp464\"",
        sep = "\n"
    )
}

bg_zenodo_metadata <- function(version, publication_date = Sys.Date()) {
    list(metadata = list(
        upload_type = "dataset",
        publication_date = format(as.Date(publication_date), "%Y-%m-%d"),
        title = sprintf("PscanR precomputed promoter backgrounds, version %d", version),
        creators = list(
            list(
                name = "Zambelli, Federico",
                affiliation = "University of Milan",
                orcid = "0000-0003-3487-4331"
            ),
            list(
                name = "Pavesi, Giulio",
                affiliation = "University of Milan",
                orcid = "0000-0001-5705-6249"
            )
        ),
        description = paste(
            "Precomputed promoter-background score distributions for the",
            "PscanR transcription factor binding motif enrichment package.",
            "This immutable release contains 105 validated combinations of",
            "JASPAR 2020, 2022, and 2024 motifs, seven genome assemblies, and",
            "five promoter windows."
        ),
        access_right = "open",
        license = "cc-by-4.0",
        version = as.character(version),
        keywords = c(
            "PscanR", "transcription factor binding sites", "motif enrichment",
            "promoters", "JASPAR", "background distributions"
        ),
        related_identifiers = list(
            list(
                identifier = "https://doi.org/10.1093/nar/gkp464",
                relation = "isSupplementTo",
                resource_type = "publication-article"
            ),
            list(
                identifier = "https://github.com/Federico77z/PscanR_backgrounds",
                relation = "isSupplementTo",
                resource_type = "dataset"
            )
        )
    ))
}

bg_write_zenodo_json <- function(metadata, path) {
    bg_require("jsonlite")
    jsonlite::write_json(
        metadata, path, auto_unbox = TRUE, pretty = TRUE, null = "null"
    )
    invisible(path)
}

bg_write_manifest <- function(root, members, path) {
    hashes <- vapply(file.path(root, members), bg_hash_file, character(1))
    writeLines(sprintf("%s  %s", hashes, members), path)
    invisible(path)
}

bg_parse_manifest <- function(path) {
    lines <- readLines(path, warn = FALSE)
    pattern <- "^([0-9a-f]{64})  (.+)$"
    if (!length(lines) || any(!grepl(pattern, lines))) {
        bg_stop("Invalid SHA-256 manifest: ", path)
    }
    data.frame(
        sha256 = sub(pattern, "\\1", lines),
        path = sub(pattern, "\\2", lines),
        stringsAsFactors = FALSE
    )
}

bg_validate_zenodo_release <- function(archive, version) {
    if (!file.exists(archive)) bg_stop("Missing Zenodo archive: ", archive)
    extract_root <- tempfile("pscanr-zenodo-validate-")
    dir.create(extract_root)
    on.exit(unlink(extract_root, recursive = TRUE), add = TRUE)
    members <- utils::unzip(archive, list = TRUE)$Name
    prefix <- sprintf("PscanR_backgrounds_v%d/", version)
    if (!length(members) || any(grepl("(^|/)\\.\\.?(/|$)", members)) ||
        any(!startsWith(members, prefix))) {
        bg_stop("Zenodo archive contains unsafe or unexpected member paths")
    }
    utils::unzip(archive, exdir = extract_root)
    release_root <- file.path(extract_root, sub("/$", "", prefix))
    catalog_path <- file.path(release_root, "catalog.tsv")
    manifest_path <- file.path(release_root, "MANIFEST.sha256")
    if (!file.exists(catalog_path) || !file.exists(manifest_path)) {
        bg_stop("Zenodo archive is missing catalog.tsv or MANIFEST.sha256")
    }
    catalog <- utils::read.delim(
        catalog_path, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character()
    )
    missing <- setdiff(BG_CATALOG_COLUMNS, names(catalog))
    if (length(missing)) {
        bg_stop("Zenodo catalog is missing columns: ", paste(missing, collapse = ", "))
    }
    if (any(catalog$status != "validated") ||
        any(catalog$background_version != as.character(version)) ||
        anyDuplicated(catalog$artifact)) {
        bg_stop("Zenodo catalog contains invalid release entries")
    }
    manifest <- bg_parse_manifest(manifest_path)
    if (anyDuplicated(manifest$path) || any(grepl("^/|(^|/)\\.\\.?(/|$)", manifest$path))) {
        bg_stop("Zenodo manifest contains duplicate or unsafe paths")
    }
    expected <- sort(setdiff(
        members[!endsWith(members, "/")],
        paste0(prefix, "MANIFEST.sha256")
    ))
    expected <- sub(paste0("^", prefix), "", expected)
    if (!identical(sort(manifest$path), expected)) {
        bg_stop("Zenodo manifest does not cover every archive payload file")
    }
    actual <- vapply(
        file.path(release_root, manifest$path), bg_hash_file, character(1)
    )
    if (!identical(unname(actual), manifest$sha256)) {
        bg_stop("Zenodo archive member failed SHA-256 validation")
    }
    artifact_paths <- file.path(release_root, catalog$artifact)
    artifact_hashes <- vapply(artifact_paths, bg_hash_file, character(1))
    if (!identical(unname(artifact_hashes), catalog$artifact_sha256)) {
        bg_stop("Zenodo catalog artifact checksum mismatch")
    }
    invisible(list(
        archive = normalizePath(archive),
        archive_sha256 = bg_hash_file(archive),
        resource_count = nrow(catalog),
        member_count = length(expected) + 1L
    ))
}

bg_build_zenodo_release <- function(root, catalog, version = 2L,
                                    output_dir = file.path(
                                        root, "releases", "zenodo",
                                        paste0("v", version)
                                    ), publication_date = Sys.Date()) {
    version <- suppressWarnings(as.integer(version))
    if (length(version) != 1L || is.na(version) || version < 1L) {
        bg_stop("Zenodo release version must be one positive integer")
    }
    release_catalog <- bg_zenodo_release_catalog(catalog, version)
    config <- bg_read_config(root)
    expected_keys <- nrow(config$organisms) * nrow(config$windows) *
        nrow(config$jaspar)
    if (nrow(release_catalog) != expected_keys) {
        bg_stop(
            "Zenodo version ", version, " has ", nrow(release_catalog),
            " resources; expected ", expected_keys
        )
    }
    source_paths <- file.path(root, "BG_files", basename(release_catalog$artifact))
    missing <- source_paths[!file.exists(source_paths)]
    if (length(missing)) {
        bg_stop("Missing release artifacts: ", paste(basename(missing), collapse = ", "))
    }
    disk_artifacts <- list.files(
        file.path(root, "BG_files"),
        pattern = sprintf("\\.psbg%d\\.txt$", version),
        full.names = TRUE
    )
    extra <- setdiff(basename(disk_artifacts), basename(source_paths))
    if (length(extra)) {
        bg_stop(
            "Unregistered release artifacts: ", paste(sort(extra), collapse = ", ")
        )
    }
    actual_hashes <- vapply(source_paths, bg_hash_file, character(1))
    if (!identical(unname(actual_hashes), release_catalog$artifact_sha256)) {
        bg_stop("One or more source artifacts do not match catalog SHA-256 values")
    }

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    archive_name <- sprintf("PscanR_backgrounds_v%d.zip", version)
    archive <- file.path(output_dir, archive_name)
    checksum_path <- paste0(archive, ".sha256")
    metadata_path <- file.path(output_dir, "zenodo_metadata.json")
    instructions_path <- file.path(output_dir, "UPLOAD_INSTRUCTIONS.md")
    unlink(c(archive, checksum_path, metadata_path, instructions_path))

    staging_parent <- tempfile("pscanr-zenodo-build-")
    release_name <- sprintf("PscanR_backgrounds_v%d", version)
    release_root <- file.path(staging_parent, release_name)
    backgrounds <- file.path(release_root, "backgrounds")
    dir.create(backgrounds, recursive = TRUE)
    on.exit(unlink(staging_parent, recursive = TRUE), add = TRUE)
    copied <- file.copy(source_paths, backgrounds, overwrite = FALSE)
    if (!all(copied)) bg_stop("Failed to stage one or more background files")
    utils::write.table(
        release_catalog, file.path(release_root, "catalog.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE, na = ""
    )
    writeLines(bg_zenodo_readme(version, release_catalog), file.path(release_root, "README.md"))
    writeLines(bg_zenodo_license(), file.path(release_root, "LICENSE"))
    writeLines(bg_zenodo_citation(version), file.path(release_root, "CITATION.cff"))
    payload <- sort(c(
        "catalog.tsv", "README.md", "LICENSE", "CITATION.cff",
        file.path("backgrounds", basename(source_paths))
    ))
    bg_write_manifest(release_root, payload, file.path(release_root, "MANIFEST.sha256"))
    all_files <- list.files(release_root, recursive = TRUE, all.files = TRUE,
                            no.. = TRUE, full.names = TRUE)
    Sys.setFileTime(all_files, BG_ZENODO_ARCHIVE_DATE)
    members <- file.path(release_name, sort(c(payload, "MANIFEST.sha256")))
    old <- setwd(staging_parent)
    on.exit(setwd(old), add = TRUE)
    status <- utils::zip(archive, files = members, flags = "-X -q")
    setwd(old)
    if (!identical(status, 0L) || !file.exists(archive)) {
        bg_stop("Failed to create Zenodo ZIP archive")
    }
    archive_hash <- bg_hash_file(archive)
    writeLines(sprintf("%s  %s", archive_hash, archive_name), checksum_path)
    bg_write_zenodo_json(
        bg_zenodo_metadata(version, publication_date), metadata_path
    )
    writeLines(c(
        sprintf("# Upload PscanR backgrounds version %d to Zenodo", version),
        "",
        sprintf("1. Upload `%s` as the only data file in a new Zenodo dataset record.", archive_name),
        "2. Enter the fields from `zenodo_metadata.json` in the Zenodo form.",
        "3. Keep file access public and the license set to CC BY 4.0.",
        "4. Before publishing, compare the uploaded file checksum with:",
        "",
        "   ```sh",
        sprintf("   sha256sum %s", archive_name),
        "   ```",
        "",
        sprintf("   Expected SHA-256: `%s`", archive_hash),
        "",
        "5. Publish the record, then record its version-specific record ID and DOI.",
        "6. Do not use the concept DOI as the ExperimentHub download target."
    ), instructions_path)
    validation <- bg_validate_zenodo_release(archive, version)
    bg_message(
        "Built %s with %d backgrounds (%s)", archive,
        validation$resource_count, validation$archive_sha256
    )
    invisible(c(validation, list(
        checksum = checksum_path,
        metadata = metadata_path,
        instructions = instructions_path
    )))
}
