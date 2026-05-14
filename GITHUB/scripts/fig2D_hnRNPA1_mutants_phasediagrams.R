## fig2D_hnRNPA1_mutants_phasediagrams.R
## Reproduces Fig 2D: 4 phase-diagram scatters of hnRNPA1 aromatic mutants
## (Aro+, WT, Aro-, Aro--) in micDROP-green.
## Source: medmax_plot.R, Loop 4 (mut.list.full) on the Fig2/2D plate.
##
## Iteration is driven by an explicit list of well IDs in the desired
## panel order. The panel title is read from the platedef row for that
## well — so the platedef CSV is the single source of truth for well-
## to-IDR mapping.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("Fig2", "2D"))
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order (matches published Fig 2D: Aro+, WT, Aro-, Aro--).
wells <- c("E10", "E01", "E11", "E13")

# Required globals for medmax.plot()
y.max <- c(); y.min <- c(); LLPS.reg.list <- list()
sat.qual <- c(); sat.val <- c(); frac.punc <- c()
plot.list <- list()
colors <- c("lightgreen", "darkgreen")  # GFP-side palette

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in Fig2/2D raw_list. Skipping.")
    next
  }
  G <- platedef$GFP[platedef$WELL == W][1]   # label from platedef
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

panel <- ggarrange(plotlist = plot.list, ncol = 1, nrow = length(wells))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(out.dir, "fig2D_hnRNPA1_mutants_phasediagrams.pdf"),
    width = 4, height = 4 * length(wells))
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "fig2D_hnRNPA1_mutants_phasediagrams.pdf"), "\n")
