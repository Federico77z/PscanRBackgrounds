# PscanRBackgrounds

Precomputed promoter backgrounds for the
[PscanR](https://github.com/Federico77z/PscanR) motif enrichment package,
distributed through Bioconductor's ExperimentHub.

PscanR compares the motif scores of a set of foreground promoters against the
distribution of the same scores over all promoters of an organism. This
package provides that distribution for 105 combinations of:

- JASPAR CORE releases 2020, 2022 and 2024;
- genome assemblies hg38, hs1, mm10, mm39, dm6, sacCer3 and TAIR9;
- promoter windows `200u_50d`, `450u_50d`, `500u_0d`, `950u_50d` and
  `1000u_0d` (bp upstream and downstream of the TSS).

Each background is a separate ExperimentHub resource. The files are hosted on
Zenodo, in one record per JASPAR release:
[2020](https://doi.org/10.5281/zenodo.23183695),
[2022](https://doi.org/10.5281/zenodo.23183713) and
[2024](https://doi.org/10.5281/zenodo.23183720). All 105 files are also
archived as a single ZIP in <https://doi.org/10.5281/zenodo.21821764>.

## Installation

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("PscanRBackgrounds")
```

## Usage

```r
library(ExperimentHub)
eh <- ExperimentHub()
query(eh, c("PscanRBackgrounds", "hg38", "JASPAR2020"))
```

PscanR retrieves and verifies the background it needs by itself; see the
PscanR vignettes.

## Generation

The pipeline that generates the backgrounds, with its configuration,
annotation snapshots and validation reports, is maintained separately at
[PscanRBackgrounds-pipeline](https://github.com/Federico77z/PscanRBackgrounds-pipeline)
and is not part of this package.

## License

The backgrounds are distributed under CC BY 4.0. Please cite the Zenodo record
and Zambelli F, Pesole G, Pavesi G (2009) Pscan: finding over-represented
transcription factor binding site motifs in sequences from co-regulated or
co-expressed genes. Nucleic Acids Research, doi:10.1093/nar/gkp464.
