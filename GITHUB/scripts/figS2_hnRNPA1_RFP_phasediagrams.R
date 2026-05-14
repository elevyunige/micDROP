## figS2_hnRNPA1_RFP_phasediagrams.R
## Reproduces Fig S2: hnRNPA1 aromatic mutants (Aro+, WT, Aro-) in micDROP-red.
## Source: medmax_plot.R, Loop 5 (mut.list.full, RFP-side) on the FigS2 plate.
##
## Note: the FigS2 platedef contains only 3 wells (Aro+, WT, Aro-);
## the Aro-- variant is not present in this plate.
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS2")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order (matches published Fig S2: Aro+, WT, Aro-).
# Per FigS2 platedef: I10=Aro+, I01=WT, I13=Aro-.
wells <- c("I10", "I01", "I13")

y.max <- c(); y.min <- c(); LLPS.reg.list <- list()
sat.qual <- c(); sat.val <- c(); frac.punc <- c()
plot.list <- list()
colors <- c("salmon1", "firebrick4")  # RFP-side palette

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS2 raw_list. Skipping.")
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

panel <- ggarrange(plotlist = plot.list, ncol = length(wells), nrow = 1)

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(out.dir, "figS2_hnRNPA1_RFP_phasediagrams.pdf"),
    width = 4 * length(wells), height = 4)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS2_hnRNPA1_RFP_phasediagrams.pdf"), "\n")
