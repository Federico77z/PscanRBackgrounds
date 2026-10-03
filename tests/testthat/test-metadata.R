metadata_file <- function() {
    system.file("extdata", "metadata.csv", package = "PscanRBackgrounds")
}
catalog_file <- function() {
    system.file("extdata", "catalog_v2.tsv", package = "PscanRBackgrounds")
}

test_that("metadata has one valid record per catalogued background", {
    metadata <- utils::read.csv(metadata_file(), stringsAsFactors = FALSE)
    catalog <- utils::read.delim(catalog_file(), stringsAsFactors = FALSE)
    expect_identical(nrow(metadata), 105L)
    expect_identical(nrow(catalog), 105L)
    expect_false(anyDuplicated(metadata$Title) > 0L)
    expect_false(anyDuplicated(metadata$RDataPath) > 0L)
    expect_setequal(
        basename(metadata$RDataPath), basename(catalog$artifact)
    )
    expect_true(all(startsWith(metadata$RDataPath, "PscanRBackgrounds/v2/")))
    expect_true(all(metadata$DispatchClass == "FilePath"))
    expect_true(all(metadata$SourceType == "TXT"))
    expect_false(anyNA(metadata$Species))
    expect_false(anyNA(metadata$TaxonomyId))
    expect_setequal(
        unique(metadata$Genome),
        c("hg38", "hs1", "mm10", "mm39", "dm6", "sacCer3", "TAIR9")
    )
})

test_that("species and taxonomy identifiers agree with the assembly", {
    metadata <- utils::read.csv(metadata_file(), stringsAsFactors = FALSE)
    expected <- c(
        hg38 = 9606L, hs1 = 9606L, mm10 = 10090L, mm39 = 10090L,
        dm6 = 7227L, sacCer3 = 4932L, TAIR9 = 3702L
    )
    expect_identical(
        unname(expected[metadata$Genome]), as.integer(metadata$TaxonomyId)
    )
})

test_that("titles encode the catalogued background", {
    metadata <- utils::read.csv(metadata_file(), stringsAsFactors = FALSE)
    stem <- sub("\\.psbg2\\.txt$", "", basename(metadata$RDataPath))
    expect_identical(metadata$Title, paste0("PscanR_bg_v2_", stem))
})
