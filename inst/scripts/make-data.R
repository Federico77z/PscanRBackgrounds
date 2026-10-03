# PscanR background resource v2: reproducible preparation recipe.
# Assisted-by: OpenAI Codex (documentation of the existing pipeline).
# This file intentionally describes the pipeline without running a genome-wide
# scan when sourced. See DATA_SOURCES.md in the installed inst directory.
#
# Immutable resource: https://doi.org/10.5281/zenodo.21821764
# Generation code and annotation snapshots:
# https://github.com/Federico77z/PscanRBackgrounds/tree/7516eee
# Scanner commit recorded in every v2 catalog row:
# https://github.com/Federico77z/PscanR/tree/47b59b6a3135a032e02fbd4eb5608397c5cd9005
#
# Inputs are specified in config/backgrounds.tsv, jaspar.tsv and windows.tsv.
# 1. Load the pinned UCSC ncbiRefSeqCurated annotations for hg38/mm10/mm39/dm6;
#    NCBI RefSeq GCF_009914755.1 for hs1; TAIR9 GFF3 for Arabidopsis; and
#    UCSC sgdGene for sacCer3. The repository annotation_snapshots directory
#    preserves transcript identifiers, chromosome, strand and TSS coordinates.
# 2. Extract strand-aware promoters from the matching BSgenome assemblies,
#    using windows 200/50, 450/50, 500/0, 950/50 and 1000/0 (upstream/downstream).
#    Restrict to configured canonical chromosomes, exclude out-of-bounds
#    promoters, retain the intended width and remove sequences >50% N.
#    Collapse identical promoter sequences, retaining the transcript legend.
# 3. Obtain JASPAR2020/2022/2024 CORE PFMatrixList objects for the matching
#    taxonomic group (vertebrates, insects, plants or fungi).
# 4. ps_build_bg() retains the best normalized PWM window over both strands
#    for each unique promoter. Ties keep the first window on a strand and
#    prefer the forward strand. The short table records motif identifier,
#    promoter count, score mean and sample standard deviation.
# 5. Validate schemas, motif IDs, ranges and SHA-256 checksums. Compare new
#    backgrounds to the previous immutable release before marking validated.
# 6. Record annotation, promoter and motif hashes, source retrieval dates,
#    scoring specification, scanner commit and runtime versions in catalog.tsv.
# 7. Package the 105 v2 tables, normalized catalog, MANIFEST.sha256, README,
#    license and citation into a deterministic ZIP. No full PSMatrixList scans
#    or raw genome files are distributed in that ZIP.
#
# In a separate reproduction checkout with required genomes already installed:
#   Rscript scripts/backgrounds.R audit
#   Rscript scripts/backgrounds.R plan --cores=1
#   Rscript scripts/backgrounds.R all --cores=1
#   Rscript scripts/backgrounds.R zenodo --release-version=2
#
# Live annotations may have changed; exact reproduction requires the archived
# inputs and scanner above. Never overwrite an existing published release or
# execute BG_scripts (historical scripts). `zenodo` packages existing validated
# tables; it does not itself reproduce their scans.
#
# Expected archive SHA-256:
# 668b80839f2b81c7f48d11d0640a06fd58774c6cc9702ec73f4f54ebe9d0b625
