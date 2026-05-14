## figS8b_6IDRs_valency_boxplots.R
## Reproduces Fig S8 panel b: 6 boxplots of Median GFP intensity by quintile
## of Valency-RFP expression, with p(%) and n labels above each panel.
##
## This is the GFP-side mirror of fig3D_10IDRs_valency_boxplots.R.
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS8")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order (matches published Fig S8b:
# FUS, hnRNPA1, HspB8, RBM14, TAF15, TDP43).
wells <- c("F23", "I23", "N23", "B23", "H23", "M23")

quintile.colors <- c("lightgoldenrod1","orange1","firebrick1",
                     "darkred","darkmagenta")
plot.list <- list()

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS8 raw_list. Skipping.")
    next
  }
  G <- platedef$GFP[platedef$WELL == W][1]
  dtmp <- data.KPlib[[W]]
  if (nrow(dtmp) == 0) next

  dtmp$bin.var <- as.numeric(dtmp$RFP_int_mean)
  brks <- quantile(dtmp$bin.var, probs = seq(0, 1, 0.20), na.rm = TRUE)
  if (any(duplicated(brks))) {
    warning("Could not compute 5 distinct quantile breaks for well ", W, "; skipping.")
    next
  }
  dtmp$Q <- cut(dtmp$bin.var, breaks = brks,
                labels = c("0.2","0.4","0.6","0.8","1.0"),
                include.lowest = TRUE)

  punc.GFP <- (dtmp$inGFPnfoci > 0 & dtmp$f1_inGFP_toGFPmean > 150)
  total_per_q <- table(dtmp$Q)
  punc_per_q  <- table(dtmp$Q[punc.GFP])
  qlev <- levels(dtmp$Q)
  p_pct <- round(as.numeric(punc_per_q[qlev]) /
                 as.numeric(total_per_q[qlev]) * 100, 1)
  n_q   <- as.numeric(punc_per_q[qlev])

  dtmp_p <- dtmp[punc.GFP, ]
  dtmp_p$GFP_int_b5 <- as.numeric(dtmp_p$GFP_int_b5)

  label_df <- data.frame(
    Q = factor(qlev, levels = qlev),
    p = ifelse(is.na(p_pct), "NA", as.character(p_pct)),
    n = ifelse(is.na(n_q),   "0",  as.character(n_q))
  )

  p <- ggplot(subset(dtmp_p, !is.na(Q)), aes(x = Q, y = GFP_int_b5, fill = Q)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Q), size = 1.5, width = 0.1, alpha = 0.6) +
    geom_text(data = label_df, aes(x = Q, y = 4800, label = p),
              inherit.aes = FALSE, size = 3) +
    geom_text(data = label_df, aes(x = Q, y = 4500, label = n),
              inherit.aes = FALSE, size = 3) +
    annotate("text", x = 0.4, y = 4800, label = "p (%):", size = 3, hjust = 0) +
    annotate("text", x = 0.4, y = 4500, label = "n =",    size = 3, hjust = 0) +
    scale_fill_manual(values = quintile.colors) +
    scale_color_manual(values = quintile.colors) +
    scale_y_continuous(limits = c(0, 5000)) +
    labs(title = G, x = "Quintile", y = "Median GFP Intensity (arb. units)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, hjust = 0.5))

  plot.list[[length(plot.list) + 1]] <- p
}

panel <- ggarrange(plotlist = plot.list, ncol = 3, nrow = 2)

pdf(file.path(out.dir, "figS8b_6IDRs_valency_boxplots.pdf"),
    width = 12, height = 7)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS8b_6IDRs_valency_boxplots.pdf"), "\n")
