# micDROP manuscript figures

R code reproducing the quantitative figure panels of:

> **Gilat, Alexandrov, Dubreuil, Levy.** Mapping interactions between disordered 
regions reveals promiscuity in biomolecular condensate formation

This repository contains scripts that re-generate figure panels from the paper.

## Repository layout

```
micDROP-figures/
├── R/
│   ├── 00_setup.R                  loads libraries + helpers, defines data.root
│   ├── 01_run.R                    sources every runner in figure order
│   ├── functions_repository.R      core plotting / Csat estimation functions
│   └── MicroscopeToolBox_minimal.R one helper kept from the original toolbox
├── scripts/
│   ├── fig2D_…  fig2E_…  fig2F_…
│   ├── fig3C_…  fig3D_…  fig3EF_…
│   ├── fig4D_…  fig4F_…
│   └── figS2_… figS3b_… figS4_… figS5_… figS6_… figS8a_… figS8b_…
│       figS9c_…  figS11A_…
└── output/                         PDFs are written here
```

## Quick start

1. **Get the data.** Per-figure raw data (one folder per figure with a
   `platedef_idr.csv` and a `results/raw_list.RDS`) is hosted separately —
   see the data-availability statement of the published paper for the
   download link. Place the entire `DATA/` tree somewhere on disk.

2. **Point the code at your data.** Open `R/00_setup.R` and edit:
   ```r
   data.root <- "/path/to/the/DATA"
   ```
   The expected layout under `data.root/` is documented at the top of that
   file.

3. **Install R packages:** `ggplot2`, `cowplot`, `ggpubr`, `tidyverse`,
   `dplyr`, `pheatmap`, `gridExtra`, `grid`, `reshape2`, `RColorBrewer`.

4. **Regenerate every figure at once.** From an R session at the
   repository root:
   ```r
   source("R/01_run.R")
   ```
   This sources every runner in `scripts/` in canonical figure order
   (Fig 2 D → 2 E → … → S11 A), captures errors per runner so a single
   failure doesn't halt the rest, and prints a summary table at the end
   listing each runner's status (`OK` / `FAILED`) and elapsed time. All
   PDFs land in `output/`.

   Or, **run a single figure** by sourcing one runner directly:
   ```r
   source("scripts/fig2D_hnRNPA1_mutants_phasediagrams.R")
   ```

## Figure-to-script mapping

| Panel              | Script                                          |
|--------------------|-------------------------------------------------|
| Fig 2 d            | `fig2D_hnRNPA1_mutants_phasediagrams.R`         | 
| Fig 2 e            | `fig2E_4IDRs_phasediagrams.R`                   | 
| Fig 2 f            | `fig2F_csat_vs_stickiness.R`                    | 
| Fig 3 c            | `fig3C_TAF15_valency_scatter.R`                 | 
| Fig 3 d            | `fig3D_10IDRs_valency_boxplots.R`               | 
| Fig 3 e, f         | `fig3EF_csat_correlations.R`                    | 
| Fig 4 d            | `fig4D_colocalization_heatmap.R`                |
| Fig 4 f            | `fig4F_TDP43_specificity_heatmap.R`             |
| Fig S2             | `figS2_hnRNPA1_RFP_phasediagrams.R`             |
| Fig S3 b           | `figS3b_csat_green_vs_red.R`                    |
| Fig S4             | `figS4_monomeric_IDRs.R`                        |
| Fig S5             | `figS5_17IDRs_green_phasediagrams.R`            | 
| Fig S6             | `figS6_20IDRs_red_phasediagrams.R`              |
| Fig S8 a           | `figS8a_6IDRs_valency_scatter.R`                | 
| Fig S8 b           | `figS8b_6IDRs_valency_boxplots.R`               | 
| Fig S9 c           | `figS9c_TDP43_phospho_boxplots.R`               | 
| Fig S11 a          | `figS11A_TDP43_deletions_phasediagrams.R`       |
-----------------------------------------------------------------------|
