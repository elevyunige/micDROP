## fig3D_10IDRs_valency_boxplots.R
## Reproduces Fig 3D: 10 boxplot panels of Median RFP intensity by quintile
## of Valency-GFP expression. Each panel labelled with p(%) per quintile.
##
## Source: closely follows valency_dependence_boxplots.R Loop 2 (RFP-side,
## green colour palette), with two modifications relative to the original:
##   1. The IDR list is the 10 IDRs from the *top* of the script (Ddx4,
##      DYRK3, ER_a, FUS, hnRNPA1, HspB8, RBM14, TAF15, TDP43, UBQ2),
##      not the NPM1/Nephrin/HP1a/NUP98 override that Loop 2 hard-coded.
##   2. The geom_text label is "p(%): X" (computed from frac.punc per
##      quintile *before* filtering to PUNC cells), rather than the
##      original "Median: X / Count: X".
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("Fig3")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order, A01..A10 in plate order (which is
# alphabetical IDR order per the platedef and matches published Fig 3D).
wells <- sprintf("A%02d", 1:10)

quintile.colors <- c("lightgoldenrod1", "darkolivegreen2", "springgreen1",
                     "springgreen3", "darkgreen")

plot.list <- list()

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in Fig3 raw_list. Skipping.")
    next
  }
  R <- platedef$RFP[platedef$WELL == W][1]
  dtmp <- data.KPlib[[W]]
  ncell <- nrow(dtmp)
  if (ncell == 0) next

  # Bin by GFP (Valency) intensity into quintiles
  dtmp$bin.var <- as.numeric(dtmp$GFP_int_mean)
  brks <- quantile(dtmp$bin.var, probs = seq(0, 1, 0.20), na.rm = TRUE)
  if (any(duplicated(brks))) {
    warning("Could not compute 5 distinct quantile breaks for well ", W, "; skipping.")
    next
  }
  dtmp$Q <- cut(dtmp$bin.var, breaks = brks,
                labels = c("0.2","0.4","0.6","0.8","1.0"),
                include.lowest = TRUE)

  punc.RFP <- (dtmp$inRFPnfoci > 0 & dtmp$f1_inRFP_toRFPmean > 150)
  total_per_q <- table(dtmp$Q)
  punc_per_q  <- table(dtmp$Q[punc.RFP])
  p_pct <- round(as.numeric(punc_per_q[levels(dtmp$Q)]) /
                 as.numeric(total_per_q[levels(dtmp$Q)]) * 100, 1)
  p_pct[is.nan(p_pct) | is.na(p_pct)] <- NA

  dtmp_p <- dtmp[punc.RFP, ]
  dtmp_p$RFP_int_b5 <- as.numeric(dtmp_p$RFP_int_b5)

  label_df <- data.frame(
    Q = factor(c("0.2","0.4","0.6","0.8","1.0"),
               levels = c("0.2","0.4","0.6","0.8","1.0")),
    p = ifelse(is.na(p_pct), "NA", as.character(p_pct))
  )

  p <- ggplot(subset(dtmp_p, !is.na(Q)), aes(x = Q, y = RFP_int_b5, fill = Q)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Q), size = 1.5, width = 0.1, alpha = 0.6) +
    geom_text(data = label_df, aes(x = Q, y = 4800, label = p),
              inherit.aes = FALSE, size = 3) +
    annotate("text", x = 0.4, y = 4800, label = "p (%):", size = 3, hjust = 0) +
    scale_fill_manual(values = quintile.colors) +
    scale_color_manual(values = quintile.colors) +
    scale_y_continuous(limits = c(0, 5000)) +
    labs(title = R, x = "Quintile", y = "Median RFP Intensity (arb. units)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, hjust = 0.5))

  plot.list[[length(plot.list) + 1]] <- p
}

panel <- ggarrange(plotlist = plot.list, ncol = 5, nrow = 2)

pdf(file.path(out.dir, "fig3D_10IDRs_valency_boxplots.pdf"),
    width = 16, height = 7)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "fig3D_10IDRs_valency_boxplots.pdf"), "\n")
