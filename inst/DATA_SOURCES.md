# Resource license and source attribution

The version-2 PscanR background archive is distributed under CC BY 4.0:
https://creativecommons.org/licenses/by/4.0/.
Credit Federico Zambelli and Giulio Pavesi and cite the immutable
record https://doi.org/10.5281/zenodo.21821764. It contains computed promoter
score distributions and a provenance catalog, not raw genome sequences.

Motifs originate from JASPAR CORE 2020, 2022 and 2024. Cite Fornes et al.
(doi:10.1093/nar/gkz1001), Castro-Mondragon et al.
(doi:10.1093/nar/gkab1113), and Rauluseviciute et al.
(doi:10.1093/nar/gkad1059), respectively. JASPAR identifies CC BY 4.0 on
https://jaspar.elixir.no/faq/; retain collection and motif identifiers.

Promoters use UCSC/NCBI RefSeq annotations for hg38, hs1, mm10, mm39 and dm6,
TAIR9 for Arabidopsis and UCSC/SGD for sacCer3, with the matching BSgenome
assemblies. The release catalog records source versions and hashes. Source
providers retain their own attribution and terms; the archive license does
not relicense the raw source data. Before submission, maintainers must confirm
reuse terms for the archived inputs, including historical TAIR9 GFF3 material
retained in the generation repository.

See `scripts/make-data.R` for the generation recipe and immutable code links.
