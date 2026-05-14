## figS5_17IDRs_green_phasediagrams.R
## Reproduces 17 of the 20 panels in Fig S5: phase-diagram scatters of
## IDRs in micDROP-green. The plate FigS5/platedef_idr.csv contains 17
## constructs (A01-A17). The published Fig S5 shows 20 panels - the
## additional Ddx4, NUP98, UBQ2 panels were assembled in Illustrator
## from a different plate.
##
## Iteration is driven by an explicit well list (A01..A17 in plate
## order). Panel title is read from the platedef row for each well.
## Note: per the FigS5 platedef, A02 is "Ddx3" (no trailing x), A16
## is "ER_a" (with underscore), and A17 is "Crowder" (= empty Venus
## + 2CLB scaffold, labelled "Dps from S. solfataricus" in the
## published figure).

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS5")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells A01..A17 in plate order
wells <- sprintf("A%02d", 1:17)

y.max <- c(); y.min <- c(); LLPS.reg.list <- list()
sat.qual <- c(); sat.val <- c(); frac.punc <- c()
plot.list <- list()
colors <- c("lightgreen", "darkgreen")

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS5 raw_list. Skipping.")
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

panel <- ggarrange(plotlist = plot.list, ncol = 4, nrow = 5)

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(out.dir, "figS5_17IDRs_green_phasediagrams.pdf"),
    width = 16, height = 20)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS5_17IDRs_green_phasediagrams.pdf"), "\n")
cat("NOTE: only 17 of the 20 published Fig S5 panels are produced here.\n",
    "      Ddx4, NUP98 and UBQ2 panels need a different plate.\n")
