# Dependencies

## R packages

- **FactoMineR**, **PTAk**, **matrixStats** -- principal component and
  multivariate data analysis.
- **DESeq2**, **edgeR**, **MLSeq** -- Bioconductor packages for
  differential expression / classification of sequencing data (install
  via `BiocManager::install()`, not CRAN).
- **ggplot2**, **gplots**, **pheatmap**, **RColorBrewer** -- plotting
  and heatmaps.
- **keras** -- deep learning (requires a working Python/TensorFlow
  backend, as with any R `keras` installation).
- **readxl** -- reading `.xlsx` files.
- **rlang**, **plyr**, **seqinr**, **maptools**, **clusterCrit**,
  **clv** -- supporting utilities and clustering evaluation.

## Matlab

`histoneAcetylation1.m` uses `imread`, which requires Matlab's **Image
Processing Toolbox**. `est_rel_entro_HJW.m` is the third-party HJW
KL-divergence estimator -- see LICENSE. Other `.m` scripts (`xlsread`,
plotting) use only core Matlab.

## Missing input data

Most scripts expect input files that are not included in this
repository -- e.g. `healthy_miRmedNor_full947.txt`,
`Data_1_12_2017_2019.csv`, `UriWatchDataSet001.xlsx`,
`UriHealthy1.xlsx`, and dozens of `LungCancerPanel*.txt` files that
`mirAnalysisClassificationTest.r` reads from a hardcoded path inside a
local FactoMineR package installation directory
(`D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/`) rather than from
this repository.

## Hardcoded paths

Many `.R`/`.r` and `.m` files contain hardcoded absolute paths to the
original author's machine. Active (non-commented) instances are flagged
with a `# EDIT:` / `% EDIT:` comment directly above them -- update
these before running a script.
