## figS9c_TDP43_phospho_boxplots.R
## Reproduces Fig S9c: per-quintile boxplots of dilute-phase GFP intensity
## for TDP43 WT, TDP43 5StoD and TDP43 12StoD phospho variants. Each
## phospho variant is GFP-tagged and paired with Valency-RFP.
## Source: valency_dependence_boxplots.R Loop 1.
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS9")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order (matches published Fig S9c:
# TDP43_WT, TDP43_5StoD, TDP43_12StoD).
# Per FigS9 platedef: B15=TDP43_WT, H15=TDP43_5StoD, D15=TDP43_12StoD.
wells <- c("B15", "H15", "D15")

quintile.colors <- c("lightgoldenrod1","orange1","firebrick1",
                     "darkred","darkmagenta")
plot.list <- list()

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS9 raw_list. Skipping.")
    next
  }
  G <- platedef$GFP[platedef$WELL == W][1]
  dtmp <- data.KPlib[[W]]
  if (nrow(dtmp) == 0) next

  dtmp$bin.var <- as.numeric(dtmp$RFP_int_mean)
  brks <- quantile(dtmp$bin.var, probs = seq(0, 1, 0.20), na.rm = TRUE)
  if (any(duplicated(brks))) next
  dtmp$Q <- cut(dtmp$bin.var, breaks = brks,
                labels = c("0.2","0.4","0.6","0.8","1.0"),
                include.lowest = TRUE)

  punc.GFP <- (dtmp$inGFPnfoci > 0 & dtmp$f1_inGFP_toGFPmean > 150)
  dtmp_p <- dtmp[punc.GFP, ]
  dtmp_p$GFP_int_b5 <- as.numeric(dtmp_p$GFP_int_b5)

  summary_stats <- dtmp_p %>%
    group_by(Q) %>%
    summarise(median = round(median(GFP_int_b5), 1),
              count  = n())

  p <- ggplot(subset(dtmp_p, !is.na(Q)),
              aes(x = Q, y = GFP_int_b5, fill = Q)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Q), size = 1.5, width = 0.1, alpha = 0.6) +
    geom_text(data = summary_stats,
              aes(x = Q, y = 4500,
                  label = paste0("Median: ", median, "\nn = ", count)),
              inherit.aes = FALSE, size = 3, vjust = -0.2) +
    scale_fill_manual(values = quintile.colors) +
    scale_color_manual(values = quintile.colors) +
    scale_y_continuous(limits = c(0, 5000)) +
    labs(title = G, x = "Quintile (Valency-RFP intensity)",
         y = "Median GFP Intensity (arb. units)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, hjust = 0.5))

  plot.list[[length(plot.list) + 1]] <- p
}

panel <- ggarrange(plotlist = plot.list, ncol = 3, nrow = 1)

pdf(file.path(out.dir, "figS9c_TDP43_phospho_boxplots.pdf"),
    width = 14, height = 4)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS9c_TDP43_phospho_boxplots.pdf"), "\n")
