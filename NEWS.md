# PscanRBackgrounds 0.99.7

- No package changes; rebuild now that BiocFileCache 3.3.2 fixes the
  ExperimentHub download failure (Bioconductor/BiocFileCache#73) that
  broke the vignette builds of 0.99.5 and 0.99.6.

# PscanRBackgrounds 0.99.6

- No changes; rebuild attempt after the vignette failed to download
  EH10548 on the build system. The cause was later traced to
  BiocFileCache 3.3.1 (see 0.99.7), not to a transient failure.

# PscanRBackgrounds 0.99.5

- Changes from the Bioconductor review: the funder is listed in `Authors@R`,
  and the metadata tests locate the installed files with plain variables.

# PscanRBackgrounds 0.99.4

- The vignette retrieves a background through ExperimentHub, now that the
  105 records are available (EH10543-EH10647).

# PscanRBackgrounds 0.99.3

- Diana Betelli is credited in the data attribution and citation, as on
  the Zenodo records.

# PscanRBackgrounds 0.99.2

- The background files are hosted on Zenodo, one record per JASPAR release;
  the ExperimentHub metadata points at them (`Location_Prefix`,
  `RDataPath`) and targets Bioconductor 3.24.

# PscanRBackgrounds 0.99.1

- Documentation: AI-assistance statement and TAIR9 reuse terms.

# PscanRBackgrounds 0.99.0

- Initial Bioconductor submission: ExperimentHub metadata and documentation
  for the 105 PscanR version-2 promoter backgrounds, one resource per
  background.
