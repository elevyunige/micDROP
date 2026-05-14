## fig3EF_csat_correlations.R
## Reproduces Fig 3E (rank-rank Csat scatter, R = 0.87) and Fig 3F (Δp vs
## Csat, R = 0.81) from correlation_plots.R block 3.
##
## All input values are hard-coded — they are summary statistics derived
## from the medmax-plot fits (Csat) and from valency_dependence (per-quintile
## Q1Csat..Q5Csat and PQ1..PQ5).
##
## Fig 3E uses base R plot()+text() exactly as in the original script.
## Fig 3F uses ggscatter with reg.line + conf.int, identical to the original.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

P     <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43","Ddx4","DYRK3","ERa","UBQ2")
Csat  <- c(  780,    271,    204,    148,    306,    442,    752,   2275,  2242,    81)
Q1Csat<- c(  679,    420,    449,    118,    442,    991,    409,   1801,   856,    48)
Q5Csat<- c( 3620,   3606,   3452,    234,   2998,   1642,   2175,  14647, 24139,   122)
Q3Csat<- c( 1459,   1535,   1500,    142,   1199,   1163,    967,   5252, 10495,    70)
Q4Csat<- c( 2413,   2404,   2010,    174,   1841,   1153,  10000,  20000, 60000,    92)
PQ1   <- c( 25.7,   36.8,   33.7,   56.7,     54,   47.9,   42.9,   19.3,  36.2,  53.9)
PQ5   <- c(  2.2,    0.6,    2.7,   55.1,    4.2,   10.2,    1.7,    0.2,   0.2,    36)
PQ3   <- c(   12,    5.8,   16.5,     53,   20.4,   37.4,     14,    3.7,   2.1,  47.4)
PQ4   <- c(  4.5,      3,    6.9,   59.3,   10.5,   28.8,    4.4,    1.1,   0.7,  44.7)

df <- data.frame(P, Csat, Q1Csat, Q5Csat, Q3Csat, Q4Csat, PQ1, PQ5, PQ3, PQ4)
df$Pdiff      <- df$PQ1 / df$PQ4
df$Csatrank   <- rank(df$Csat)
df$Q4Csatrank <- rank(df$Q4Csat)

# ---- Fig 3E: rank-vs-rank, base R plot + text -----------------------
pdf(file.path(out.dir, "fig3E_csat_rank_corr.pdf"), width = 6, height = 6)
plot(df$Csatrank, df$Q4Csatrank,
     xlab = "Csat (rank)", ylab = "Q0.8 Csat (rank)",
     pch = 19, col = "darkgreen")
text(df$Csatrank, df$Q4Csatrank, lab = df$P, pos = 4, cex = 0.9)
abline(lm(Q4Csatrank ~ Csatrank, data = df), col = "darkgreen")
r <- cor.test(df$Csatrank, df$Q4Csatrank, method = "pearson")
mtext(sprintf("R = %.2f, p = %.4f", r$estimate, r$p.value),
      side = 3, line = 0.3, cex = 0.9)
dev.off()
cat("Wrote:", file.path(out.dir, "fig3E_csat_rank_corr.pdf"), "\n")

# ---- Fig 3F: ggscatter Pdiff vs Csat, log-log -----------------------
sp <- ggscatter(df, x = "Csat", y = "Pdiff", label = "P",
                add = "reg.line",
                add.params = list(color = "darkgreen", fill = "lightgray"),
                conf.int = TRUE)

p <- sp +
  stat_cor(method = "pearson", label.x = log10(5), label.y = log10(1)) +
  scale_x_continuous(trans = "log", breaks = c(1, 10, 100, 500, 1000, 2000)) +
  scale_y_continuous(trans = "log", breaks = c(1, 10, 100))

pdf(file.path(out.dir, "fig3F_pdiff_vs_csat.pdf"), width = 6, height = 5)
print(p)
dev.off()
cat("Wrote:", file.path(out.dir, "fig3F_pdiff_vs_csat.pdf"), "\n")
