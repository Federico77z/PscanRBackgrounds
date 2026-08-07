#' Precomputed promoter backgrounds for PscanR
#'
#' PscanRBackgrounds provides ExperimentHub metadata and documentation for
#' immutable promoter-background score distributions used by PscanR.
#'
#' @keywords internal
#' @name PscanRBackgrounds-package
NULL

#' PscanR promoter backgrounds, version 2
#'
#' A ZIP archive containing 105 validated promoter-background score
#' distributions for combinations of JASPAR 2020, 2022, and 2024, seven genome
#' assemblies, and five promoter windows. The archive also contains the
#' versioned resource catalog, SHA-256 manifest, citation information, and
#' format documentation.
#'
#' The resource uses ExperimentHub's `FilePath` dispatch and therefore returns
#' the path to the cached ZIP archive. PscanR selects and validates the required
#' background within the archive.
#'
#' @format A cached file path to `PscanR_backgrounds_v2.zip`.
#' @return A character scalar containing the path to the cached ZIP archive.
#' @source \doi{10.5281/zenodo.21821764}
#' @references Zambelli F, Pesole G, Pavesi G. (2009). Pscan: finding
#'   over-represented transcription factor binding site motifs in sequences
#'   from co-regulated or co-expressed genes. Nucleic Acids Research.
#'   \doi{10.1093/nar/gkp464}.
#' @examples
#' if (interactive()) {
#'     archive <- PscanR_backgrounds_v2()
#'     archive
#' }
#' @name PscanR_backgrounds_v2
NULL
