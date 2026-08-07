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

## ExperimentHub distribution

This repository is also the source of the lightweight `PscanRBackgrounds`
ExperimentHub package. The validated version-2 collection is archived at
[Zenodo](https://doi.org/10.5281/zenodo.21821764). The package metadata points
to that immutable version-specific record; GitHub files remain a legacy
fallback and are excluded from the package tarball.

Validate the Hub metadata and build the lightweight source package with:

```sh
Rscript inst/scripts/make-metadata.R
R CMD build .
R CMD check --as-cran --no-manual ../PscanRBackgrounds_0.99.0.tar.gz
```

The metadata targets Bioconductor 3.24 and dispatches the Zenodo ZIP as a
cached `FilePath`. A PDF manual check additionally requires `pdflatex`.

## Runtime dependencies

Use the same Bioconductor release as the PscanR checkout. The pipeline requires
PscanR plus `pkgload`, `txdbmaker`, `GenomicFeatures`, `GenomeInfoDb`,
`Biostrings`, `BSgenome`, `rtracklayer`, `DBI`, `RMariaDB`, `RSQLite`,
`JASPAR2020`, `JASPAR2022`, and `JASPAR2024`. Install the supported genomes
before a complete run:

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

Compare two complete immutable background versions before promoting a new
release:

```sh
Rscript scripts/backgrounds.R compare \
  --reference-version=1 \
  --candidate-version=2 \
  --benchmark-dir=../Test/pscan_benchmark/baseline_independent_current
```

The comparison checks all paired motif profiles, quantifies promoter-universe
changes, and recalculates foreground Z scores, ranks, and FDR values from the
stored benchmark tables without rescanning DNA. Reports and acceptance checks
are written to `reports/comparison_v1_v2/`; the command exits unsuccessfully
when a scientific compatibility threshold is exceeded. The benchmark argument
is optional, but should be supplied for release validation.

Build the validated archive and metadata for an immutable Zenodo release:

```sh
Rscript scripts/backgrounds.R zenodo --release-version=2
```

The command selects the complete validated release from `catalog.tsv`, checks
every source checksum, creates a deterministic ZIP, extracts it again, and
validates its catalog and manifest. Upload-ready files are written under
`releases/zenodo/v2/`. Only the ZIP is deposited as the Zenodo data file; the
generated metadata and instructions describe the accompanying record.

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

The calibration keeps set sizes 5, 10, and 20 as small-set diagnostics and
gates release compatibility from 50 promoters upward. Override the boundary with
`--minimum-set-size=N`. False-positive-rate uncertainty is estimated across
the independently sampled promoter sets, preserving correlation among motifs;
per-motif KS uniformity remains a diagnostic rather than a release gate.

## Version and regeneration rules

- Versions are immutable positive integers scoped to one background key.
- Only validated entries may be marked `latest`.
- A new version must pass the direct profile comparison and downstream
  benchmark checks before it is published as the recommended release.
- Regeneration occurs only when the final unique promoter-sequence hash, motif
  content hash, or explicit scoring specification changes.
- Annotation or package metadata changes that leave computational inputs
  unchanged do not cause scans.
- Performance-only PscanR changes do not require regeneration when scanner
  compatibility tests establish identical scores.
