test_that("ExperimentHub metadata identifies the immutable Zenodo release", {
    path <- system.file(
        "extdata", "metadata.csv", package = "PscanRBackgrounds"
    )
    expect_true(nzchar(path))

    metadata <- utils::read.csv(path, stringsAsFactors = FALSE)
    expect_identical(nrow(metadata), 1L)
    expect_identical(metadata$Title, "PscanR_backgrounds_v2")
    expect_identical(metadata$SourceVersion, 2L)
    expect_identical(
        metadata$SourceUrl,
        "https://doi.org/10.5281/zenodo.21821764"
    )
    expect_identical(metadata$DispatchClass, "FilePath")
    expect_identical(
        metadata$RDataPath,
        "PscanR_backgrounds_v2.zip/content"
    )
})

test_that("the installed package exports the documented Hub accessor", {
    path <- system.file("extdata", "metadata.csv", package="PscanRBackgrounds")
    title <- utils::read.csv(path)$Title
    expect_true(title %in% getNamespaceExports("PscanRBackgrounds"))
    accessor <- getExportedValue("PscanRBackgrounds", title)
    expect_true(is.function(accessor))
    expect_true("metadata" %in% names(formals(accessor)))
})
