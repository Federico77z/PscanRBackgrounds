# PscanR Backgrounds

This repository stores precomputed promoter-background statistics used by
[PscanR](https://github.com/Federico77z/PscanR). For every supported assembly,
promoter window, and JASPAR CORE taxonomic collection, PscanR retains the best
normalized PWM match in each unique promoter sequence and stores the resulting
mean and standard deviation.

## Supported combinations

- JASPAR releases: 2020, 2022, and 2024.
- Assemblies: hg38, hs1, mm10, mm39, dm6, sacCer3, and TAIR9.
- Promoter windows: `200u_50d`, `450u_50d`, `500u_0d`, `950u_50d`, and
  `1000u_0d`.

The complete matrix contains 105 background keys. `catalog.tsv` is the
authoritative artifact registry. Existing version-1 entries are immutable
legacy backgrounds; they predate input fingerprinting. New versions include
annotation, promoter-sequence, motif-content, scanner, and artifact hashes.

## Runtime dependencies

Use the same Bioconductor release as the PscanR checkout. The pipeline requires
PscanR plus `pkgload`, `txdbmaker`, `GenomicFeatures`, `GenomeInfoDb`,
`Biostrings`, `BSgenome`, `rtracklayer`, `RSQLite`, `JASPAR2020`, `JASPAR2022`,
and `JASPAR2024`. Install the supported genomes before a complete run:

```r
BiocManager::install(c(
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Hsapiens.UCSC.hs1",
    "BSgenome.Mmusculus.UCSC.mm10",
    "BSgenome.Mmusculus.UCSC.mm39",
    "BSgenome.Dmelanogaster.UCSC.dm6",
    "BSgenome.Scerevisiae.UCSC.sacCer3",
    "BSgenome.Athaliana.TAIR.TAIR9"
))
```

The command fails before scanning if a selected genome or package is absent.

## Master command

Run commands from the repository root:

```sh
Rscript scripts/backgrounds.R audit
Rscript scripts/backgrounds.R check
Rscript scripts/backgrounds.R plan --cores=16
Rscript scripts/backgrounds.R all --cores=16
```

`check` is a lightweight annotation-coordinate check. `plan` constructs the
exact unique promoter sequences and motif collections, calculates their hashes,
and reports which jobs are required. `all` generates only those jobs, validates
them in `staging/`, publishes them without overwriting an existing version, and
updates the catalog atomically.

Use filters to run part of the matrix:

```sh
Rscript scripts/backgrounds.R plan \
  --assembly=hg38 --jaspar=2024 --window=450u_50d
```

The sibling PscanR checkout is used automatically. Set `PSCANR_SOURCE` when it
is elsewhere. Generation refuses a dirty PscanR checkout because the exact
scanner commit is part of provenance.

## Scientific validation

The calibration is intentionally separate from routine generation because it
is expensive:

```sh
Rscript scripts/backgrounds.R calibrate --cores=16
```

It reproduces random-promoter null tests for all assemblies using JASPAR2024
and the `450u_50d` window. It also compares production backgrounds with a
deterministic one-promoter-per-gene diagnostic universe. Reports are written
under `reports/`; the command exits unsuccessfully when the configured
calibration criteria do not pass. Treat that result as a publication blocker
and review it before running `all`.

## Version and regeneration rules

- Versions are immutable positive integers scoped to one background key.
- Only validated entries may be marked `latest`.
- Regeneration occurs only when the final unique promoter-sequence hash, motif
  content hash, or explicit scoring specification changes.
- Annotation or package metadata changes that leave computational inputs
  unchanged do not cause scans.
- Performance-only PscanR changes do not require regeneration when scanner
  compatibility tests establish identical scores.
