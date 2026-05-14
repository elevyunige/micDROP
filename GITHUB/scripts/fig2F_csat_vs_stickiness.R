## fig2F_csat_vs_stickiness.R
## Reproduces Fig 2F: scatter of Csat vs Stickiness for 10 IDRs (Pearson R = -0.72).
## Source: correlation_plots.R, block 2.
##
## Csat values (RFP, micDROP-red) and Stickiness sums are hard-coded — they are
## summary statistics derived from the medmax-plot fits in Fig S6 (Csat) and from
## independent IDR sequence analysis (Stickiness, sum of pairwise residue
## interaction strengths).

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

P          <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43","Ddx4","DYRK3","ERa","UBQ2")
Csat       <- c(  780,    271,    423,    148,    306,    442,   752,   2275,  2242,    81)
Stickiness <- c(5.938,  8.964, 10.877, 20.031,  5.097,  8.576, 6.284,  8.197, 7.543, 22.408)

df <- data.frame(P, Csat, Stickiness)

sp <- ggscatter(df, x = "Stickiness", y = "Csat", label = "P",
                add = "reg.line",
                add.params = list(color = "darkgreen", fill = "lightgray"),
                conf.int = TRUE)

p <- sp +
  stat_cor(method = "pearson", label.x = 10, label.y = 10) +
  scale_y_continuous(trans = "log", breaks = c(1, 10, 100, 500, 1000, 2000))

pdf(file.path(out.dir, "fig2F_csat_vs_stickiness.pdf"), width = 6, height = 5)
print(p)
dev.off()

cat("Wrote:", file.path(out.dir, "fig2F_csat_vs_stickiness.pdf"), "\n")
