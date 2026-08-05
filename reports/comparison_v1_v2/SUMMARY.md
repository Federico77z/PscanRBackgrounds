# Background version 1 versus version 2

Overall result: **PASS**.

## Profile comparison

- Paired background keys: 105.
- Keys with changed promoter counts: 15; unchanged: 90.
- Motif additions/removals: 0/0.
- Maximum normalized mean shift in changed universes: 2.03565e-03 SD.
- Maximum projected null Z-score shift at n=1000: 6.43728e-02.
- Documented zero-variance SD sentinel changes: 13 motifs.

## Downstream benchmark

- Independent foreground datasets: 100.
- Benchmark manifest SHA-256: `9c1c706b051b8a6cd04401ca26a30f05c7c08f58716b6e19243106795ca9d9cd`.
- Maximum absolute Z-score delta: 4.81421e-02.
- Minimum rank Spearman correlation: 9.99862e-01.
- FDR classification changes: 47 of 74600 motif/dataset pairs (6.30027e-04).

## Acceptance checks

| Scope | Metric | Observed | Criterion | Passed |
|---|---|---:|---:|:---:|
| profiles | paired background keys | 105 | = 105 | yes |
| profiles | motif set differences | 0 | = 0 | yes |
| stable profiles | maximum mean delta | 1.11022302462516e-15 | <= 1e-12 | yes |
| stable profiles | maximum SD delta excluding sentinels | 9.99200722162641e-16 | <= 1e-12 | yes |
| changed profiles | minimum mean/SD correlation | 0.999997360985748 | >= 0.9999 | yes |
| changed profiles | maximum normalized mean shift | 0.00203564719254603 | <= 0.01 | yes |
| changed profiles | minimum SD ratio | 0.998125811928348 | >= 0.95 | yes |
| changed profiles | maximum SD ratio | 1.00226381558787 | <= 1.05 | yes |
| changed profiles | SD ratio first percentile | 0.999017484778889 | >= 0.98 | yes |
| changed profiles | SD ratio 99th percentile | 1.00114897264535 | <= 1.02 | yes |
| changed profiles | maximum projected null Z shift at n=1000 | 0.0643728164097279 | <= 0.25 | yes |
| benchmark | reference background value error | 0 | <= 1e-12 | yes |
| benchmark | reference Z recalculation error | 7.93143328792212e-13 | <= 1e-10 | yes |
| benchmark | minimum rank Spearman | 0.999862327840249 | >= 0.999 | yes |
| benchmark | minimum top-10 overlap | 0.9 | >= 0.8 | yes |
| benchmark | minimum top-25 overlap | 0.96 | >= 0.9 | yes |
| benchmark | minimum top-50 overlap | 0.96 | >= 0.9 | yes |
| benchmark | maximum absolute Z delta | 0.0481421273423355 | <= 0.1 | yes |
| benchmark | FDR classification discordance rate | 0.000630026809651475 | <= 0.005 | yes |
| benchmark | maximum discordant FDR boundary distance | 0.0045000010894953 | <= 0.01 | yes |
