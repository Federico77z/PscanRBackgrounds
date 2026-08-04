# Initial Background Audit

Audit date: 2026-08-04

## Confirmed behavior

The implementation follows the Pscan model: all configured promoter sequences
are scanned on both strands, the best normalized PWM match per promoter is
retained, and the promoter-universe mean and sample standard deviation are used
for foreground z-tests. The current production universe contains unique
transcript-promoter sequences, not one representative promoter per gene.

All 105 inherited tables parse successfully. Each has the expected four fields,
one positive promoter count, unique motif IDs, means in `[0, 1]`, and positive
standard deviations no greater than one. For motif IDs shared between JASPAR
2020, 2022, and 2024, background values are identical within numerical
precision for the same assembly and promoter window.

## Corrected defects

- Three valid tables used filenames that PscanR could not construct. They were
  renamed to the canonical version-1 convention.
- A machine-readable catalog and SHA-256 checksum were added for every table.
- The duplicated generation scripts were superseded by a configuration-driven,
  resumable master command.

## Legacy limitations

Version-1 generation did not preserve live UCSC annotation snapshots, source
retrieval dates, promoter hashes, motif hashes, PscanR commits, or Bioconductor
versions. Exact input-level reproduction therefore cannot be proven
retrospectively. These files remain immutable legacy artifacts. The first
fingerprinted regeneration will create version 2 where inputs have changed or
cannot be proven unchanged.

The old scripts also selected canonical chromosomes by positional index,
disabled TLS verification, and contained hs1 scripts that depended on an
interactive `gff` object. They remain in `BG_scripts/` for historical review
only and must not be used for new production backgrounds.

## Open scientific validation

The transcript-level universe is retained for compatibility. Before considering
a switch to one promoter per gene, run the master command's `calibrate` mode and
review `reports/scientific_calibration_summary.tsv` together with
`reports/promoter_universe_comparison.tsv`.

