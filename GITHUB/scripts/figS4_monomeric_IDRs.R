## figS4_monomeric_IDRs.R
## Reproduces Fig S4 panel a: 16 phase-diagram scatters of monomeric
## (non-Dps-fused) IDRs in micDROP-green.
## Source: medmax_plot.R, Loop 3 (GFPmono.list.full) on the FigS4 plate.
##
## Iteration is driven by an explicit well list (D01..D16 in plate
## order). Panel title is read from the platedef row for each well.
## Note: per the FigS4 platedef, well D03 is "GATA_mono" (no `3`)
## and D14 is "Era_mono" (no underscore). The published figure relabels
## these to GATA3 and ERa.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS4")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells D01..D16 in plate order (matches the order monomeric IDRs were
# pipetted, which is also the order shown in published Fig S4).
wells <- sprintf("D%02d", 1:16)

y.max <- c(); y.min <- c(); LLPS.reg.list <- list()
sat.qual <- c(); sat.val <- c(); frac.punc <- c()
plot.list <- list()
colors <- c("lightgreen", "darkgreen")

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS4 raw_list. Skipping.")
    next
  }
  G <- platedef$GFP[platedef$WELL == W][1]
  dtmp   <- data.KPlib[[W]]
  ncell  <- nrow(dtmp)
  strain <- paste("YFP-", G)
  punc   <- (dtmp$inGFPnfoci > 0 & dtmp$f1_inGFP_toGFPmean > 150)
  counts <- sum(punc, na.rm = TRUE)
  frac.punc[I] <- counts / ncell
  dtmp$punc <- punc

  if (ncell > 0 & counts > 50) {
    y.max[I]      <- quantile(dtmp$GFP_int_b5[punc], prob = 0.90)
    y.min[I]      <- quantile(dtmp$GFP_int_b5[punc], prob = 0.20)
    punc.GOOD     <- punc & dtmp$GFP_int_b5 < y.max[I] & dtmp$GFP_int_b5 > y.min[I]
    LLPS.reg      <- lm(log(dtmp$GFP_int_b5[punc.GOOD]) ~ log(dtmp$GFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] <- LLPS.reg
    sat.val[I]    <- mean(dtmp$GFP_int_b5[punc.GOOD])
    sat.qual[I]   <- LLPS.reg$coefficients[2]
    plot.list[[I]] <- medmax.plot(FP = "GFP", ST = 1)
  } else if (ncell > 0 & counts > 0) {
    plot.list[[I]] <- medmax.plot(FP = "GFP", ST = 2)
  } else {
    plot.list[[I]] <- medmax.plot(FP = "GFP", ST = 3)
  }
}

panel <- ggarrange(plotlist = plot.list, ncol = 4, nrow = 4)

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(out.dir, "figS4_monomeric_IDRs.pdf"), width = 16, height = 16)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS4_monomeric_IDRs.pdf"), "\n")
