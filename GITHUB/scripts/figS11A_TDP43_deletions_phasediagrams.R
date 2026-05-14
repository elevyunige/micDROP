## figS11A_TDP43_deletions_phasediagrams.R
## Reproduces Fig S11a: 4 phase-diagram scatters of TDP43 alpha-helix
## deletion mutants (WT, Δ320-330, Δ331-343, Δ320-343) in micDROP-green.
## Source: medmax_plot.R, Loop 4 (mut.list.full) on the FigS11/A plate.
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("FigS11", "A"))
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order (matches published Fig S11a:
# WT, Δ320-330, Δ331-343, Δ320-343).
# Per FigS11/A platedef: E02=TDP43_WT, E05=TDP43_320_330,
# E06=TDP43_331_343, E09=TDP43_320_343.
wells <- c("E02", "E05", "E06", "E09")

y.max <- c(); y.min <- c(); LLPS.reg.list <- list()
sat.qual <- c(); sat.val <- c(); frac.punc <- c()
plot.list <- list()
colors <- c("lightgreen", "darkgreen")

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS11/A raw_list. Skipping.")
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

panel <- ggarrange(plotlist = plot.list, ncol = length(wells), nrow = 1)

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(out.dir, "figS11A_TDP43_deletions_phasediagrams.pdf"),
    width = 4 * length(wells), height = 4)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS11A_TDP43_deletions_phasediagrams.pdf"), "\n")
