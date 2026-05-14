## figS3b_csat_green_vs_red.R
## Reproduces Fig S3 panel b: per-IDR Csat measured in micDROP-green vs
## micDROP-red, R = 0.93. Source: correlation_plots.R block 1.
##
## NOTE: Fig S3 panels a (UMAP) and c (Csat-vs-IDR-feature correlations)
## are NOT reproduced from this repo - they require additional analyses
## (UMAP with the umap library, plus IDR feature computation) that are
## not part of the scripts shipped here.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

P   <- c("DYRK3","ERa","FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43")
GFP <- c(  758,   842,  399,     53,     204,    57,    182,   477)
RFP <- c( 2275,  2242,  780,    271,     423,   148,    306,   442)
df  <- data.frame(P, GFP, RFP)

sp <- ggscatter(df, x = "GFP", y = "RFP", label = "P",
                add = "reg.line",
                add.params = list(color = "darkgreen", fill = "lightgray"),
                conf.int = TRUE)

p <- sp + stat_cor(method = "pearson", label.x = 3, label.y = 30) +
  labs(x = "micDROP-green Csat (arb. units)",
       y = "micDROP-red Csat (arb. units)")

pdf(file.path(out.dir, "figS3b_csat_green_vs_red.pdf"), width = 6, height = 5)
print(p)
dev.off()

cat("Wrote:", file.path(out.dir, "figS3b_csat_green_vs_red.pdf"), "\n")
