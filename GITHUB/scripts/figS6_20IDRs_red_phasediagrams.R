## figS6_20IDRs_red_phasediagrams.R
## Reproduces Fig S6 panel a: 20 phase-diagram scatters of IDRs in
## micDROP-red. Source: medmax_plot.R, Loop 1 (RFP.list.full).
##
## Iteration is driven by an explicit well list (A01..A20 in plate
## order). Panel title is read from the platedef row for each well.
## Note: per the FigS6 platedef, A02 is "Ddx3" (no trailing x), A18
## is "ER_a" (with underscore), and A15 is "2CLB" (the empty mScarlet
## + S. solfataricus Dps scaffold; the published figure labels this
## panel "Dps from V. cholerae").

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS6")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells A01..A20 in plate order
wells <- sprintf("A%02d", 1:20)

y.max <- c(); y.min <- c(); LLPS.reg.list <- list()
sat.qual <- c(); sat.val <- c(); frac.punc <- c()
plot.list <- list()
colors <- c("salmon1", "firebrick4")

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS6 raw_list. Skipping.")
    next
  }
  R <- platedef$RFP[platedef$WELL == W][1]
  dtmp   <- data.KPlib[[W]]
  ncell  <- nrow(dtmp)
  strain <- paste("RFP-", R)
  punc   <- (dtmp$inRFPnfoci > 0 & dtmp$f1_inRFP_toRFPmean > 150)
  counts <- sum(punc, na.rm = TRUE)
  frac.punc[I] <- counts / ncell
  dtmp$punc <- punc

  if (ncell > 0 & counts > 50) {
    y.max[I]      <- quantile(dtmp$RFP_int_b5[punc], prob = 0.90)
    y.min[I]      <- quantile(dtmp$RFP_int_b5[punc], prob = 0.20)
    punc.GOOD     <- punc & dtmp$RFP_int_b5 < y.max[I] & dtmp$RFP_int_b5 > y.min[I]
    LLPS.reg      <- lm(log(dtmp$RFP_int_b5[punc.GOOD]) ~ log(dtmp$RFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] <- LLPS.reg
    sat.val[I]    <- mean(dtmp$RFP_int_b5[punc.GOOD])
    sat.qual[I]   <- LLPS.reg$coefficients[2]
    plot.list[[I]] <- medmax.plot(FP = "RFP", ST = 1)
  } else if (ncell > 0 & counts > 0) {
    plot.list[[I]] <- medmax.plot(FP = "RFP", ST = 2)
  } else {
    plot.list[[I]] <- medmax.plot(FP = "RFP", ST = 3)
  }
}

panel <- ggarrange(plotlist = plot.list, ncol = 4, nrow = 5)

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(out.dir, "figS6_20IDRs_red_phasediagrams.pdf"),
    width = 16, height = 20)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS6_20IDRs_red_phasediagrams.pdf"), "\n")
