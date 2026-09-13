## DiseaseDetectionUrine

## Quick start

This repository contains R and Matlab scripts for companion diagnostics
(UriWatch - LiquidThermometer) based on urinary microRNA biomarkers.
See [DEPENDENCIES.md](DEPENDENCIES.md) for the required R packages and
the input data files scripts expect but which are not included in this
repository.

## Repository contents

- `mirAnalysis.r`, `mirAnalysisClassificationTest.r`,
  `mirAnalysisClassificationTest2.r`, `micro_RNA_analyses_1_12_2017.R`,
  `stepsRUriWatch.r` -- microRNA classification and dose-response
  analysis.
- `CA.R`, `PCA.R`, `Heatmap1.r`, `hclust.R`, `ColorBrewer.R`, `v25i01.R`,
  `R scripts for PONE-D-15-46623.R` -- multivariate data analysis,
  PCA, and heatmap generation.
- `est_rel_entro_HJW.m` -- the third-party HJW KL-divergence estimator
  (see LICENSE).
- `HealthyThreeTimePointMirs.m`, `HealthyThreeTimePointMirsLung.m`,
  `HealthyThreeTimePointMirsMM.m`, `uriWatchHistogram.m`,
  `histoneAcetylation1.m`, `granger_cause.m`, `miRloess.m`,
  `miRnaTrim.m`, `miRnormalize.m`, `quantile_norm.m`, `roc_curve.m` --
  Matlab analysis and utility scripts.
- `data/` -- `.RData` data files.
- [`media/`](media/) -- figures.
- `reports/` -- a PDF report.
- **License:** see [LICENSE](LICENSE) -- research/educational use.

## About

R and Matlab code I wrote for companion diagnostics, UriWatch - LiquidThermometer 

Matlab code I wrote to segregate lung cancer stage IV drug-naive patients in Rotterdam and Enschede, NL (provided by Johan de Rooij, PhD) from healthy urinary microRNA samples (the divergence implementation is by Tsachy Weissman lab) 

R code I wrote for disease biomarker multivariate data analysis (the principal component analysis is by FactoMineR, there is also a function for differential expression analysis and heatmaps from the Thomas Tuschl lab, an additional function for differential expression analysis, support vector machine and random forest classification is by Karin Groothuis-Oudshoorn, PhD) of urinary microRNA

For detailed information, see: https://www.journalofliquidbiopsy.com/article/S2950-1954(26)00002-0/fulltext

