## figS8a_6IDRs_valency_scatter.R
## Reproduces Fig S8 panel a: 6 quintile-coloured phase-diagram scatters
## of GFP-tagged IDRs (FUS, hnRNPA1, HspB8, RBM14, TAF15, TDP43) paired
## with Valency-RFP. Source: valency_dependence.R, Loop 1.
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("FigS8")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Wells in desired panel order (matches published Fig S8a:
# FUS, hnRNPA1, HspB8, RBM14, TAF15, TDP43).
# Per FigS8 platedef: F23=FUS, I23=hnRNPA1, N23=HspB8, B23=RBM14,
# H23=TAF15, M23=TDP43.
wells <- c("F23", "I23", "N23", "B23", "H23", "M23")

plot.list <- list()

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in FigS8 raw_list. Skipping.")
    next
  }
  G <- platedef$GFP[platedef$WELL == W][1]
  dtmp   <- data.KPlib[[W]]
  ncell  <- nrow(dtmp)
  strain <- paste("YFP-", G)
  dtmp$GFP_RFP_b5 <- as.numeric(dtmp$RFP_int_mean)
  if (ncell == 0) {
    plot.list[[I]] <- val.plot.full(n = 0, FP = "GFP"); next
  }

  q1 <- quantile(dtmp$GFP_RFP_b5, prob = 0.20)
  q2 <- quantile(dtmp$GFP_RFP_b5, prob = 0.40)
  q3 <- quantile(dtmp$GFP_RFP_b5, prob = 0.60)
  q4 <- quantile(dtmp$GFP_RFP_b5, prob = 0.80)
  punc.GFP <- (dtmp$inGFPnfoci > 0 & dtmp$f1_inGFP_toGFPmean > 150)

  data.q1 <- dtmp %>% filter(GFP_RFP_b5 <= q1                       & punc.GFP)
  data.q2 <- dtmp %>% filter(GFP_RFP_b5 <= q2 & GFP_RFP_b5 > q1     & punc.GFP)
  data.q3 <- dtmp %>% filter(GFP_RFP_b5 <= q3 & GFP_RFP_b5 > q2     & punc.GFP)
  data.q4 <- dtmp %>% filter(GFP_RFP_b5 <= q4 & GFP_RFP_b5 > q3     & punc.GFP)
  data.q5 <- dtmp %>% filter(GFP_RFP_b5  > q4                        & punc.GFP)
  data <- list(data.q1, data.q2, data.q3, data.q4, data.q5)

  frac.punc.q1 <- nrow(data.q1) / sum(dtmp$GFP_RFP_b5 <  q1)
  frac.punc.q2 <- nrow(data.q2) / sum(dtmp$GFP_RFP_b5 <  q2 & dtmp$GFP_RFP_b5 > q1)
  frac.punc.q3 <- nrow(data.q3) / sum(dtmp$GFP_RFP_b5 <  q3 & dtmp$GFP_RFP_b5 > q2)
  frac.punc.q4 <- nrow(data.q4) / sum(dtmp$GFP_RFP_b5 <  q4 & dtmp$GFP_RFP_b5 > q3)
  frac.punc.q5 <- nrow(data.q5) / sum(dtmp$GFP_RFP_b5 >  q4)
  frac.punc.list <- list(round(frac.punc.q1,3), round(frac.punc.q2,3),
                         round(frac.punc.q3,3), round(frac.punc.q4,3),
                         round(frac.punc.q5,3))

  ns <- c(nrow(data.q1), nrow(data.q2), nrow(data.q3), nrow(data.q4), nrow(data.q5))
  n.q <- if (all(ns > 20)) 5L else
         if (all(ns[1:4] > 20)) 4L else
         if (all(ns[1:3] > 20)) 3L else
         if (all(ns[1:2] > 20)) 2L else
         if (ns[1] > 20) 1L else 0L

  if (n.q > 0) {
    SUM <- fig.data(n = n.q, FP = "GFP")
    plot.list[[I]] <- val.plot.full(n = n.q, FP = "GFP")
  } else {
    plot.list[[I]] <- val.plot.full(n = 0, FP = "GFP")
  }
}

panel <- ggarrange(plotlist = plot.list, ncol = 3, nrow = 2)

pdf(file.path(out.dir, "figS8a_6IDRs_valency_scatter.pdf"),
    width = 18, height = 14)
print(panel)
dev.off()

cat("Wrote:", file.path(out.dir, "figS8a_6IDRs_valency_scatter.pdf"), "\n")
