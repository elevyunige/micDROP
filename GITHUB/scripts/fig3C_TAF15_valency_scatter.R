## fig3C_TAF15_valency_scatter.R
## Reproduces Fig 3C: TAF15 cells with condensate, scatter of Median RFP vs
## Maximal RFP, points colored by GFP (Valency) quintile, per-quintile
## regression lines.
## Source: valency_dependence.R, Loop 2 (RFP-tagged IDR x Valency-GFP),
## restricted to TAF15.
##
## Iteration is driven by an explicit well list. Panel title from platedef.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data("Fig3")
data.KPlib <- dat$data
platedef   <- dat$platedef

# Per Fig3 platedef, A08 is the TAF15 row (Valency-GFP / TAF15-RFP).
wells <- c("A08")

plot.list <- list()

for (I in seq_along(wells)) {
  W <- wells[I]
  if (is.null(data.KPlib[[W]])) {
    warning("No data for well '", W, "' in Fig3 raw_list. Skipping.")
    next
  }
  R <- platedef$RFP[platedef$WELL == W][1]
  dtmp   <- data.KPlib[[W]]
  ncell  <- nrow(dtmp)
  strain <- paste("RFP-", R)
  dtmp$GFP_RFP_b5 <- as.numeric(dtmp$GFP_int_mean)

  if (ncell == 0) {
    plot.list[[I]] <- val.plot.full(n = 0, FP = "RFP")
    next
  }

  q1 <- quantile(dtmp$GFP_RFP_b5, prob = 0.20)
  q2 <- quantile(dtmp$GFP_RFP_b5, prob = 0.40)
  q3 <- quantile(dtmp$GFP_RFP_b5, prob = 0.60)
  q4 <- quantile(dtmp$GFP_RFP_b5, prob = 0.80)
  punc.RFP <- (dtmp$inRFPnfoci > 0 & dtmp$f1_inRFP_toRFPmean > 150)

  data.q1 <- dtmp %>% filter(GFP_RFP_b5 <= q1                       & punc.RFP)
  data.q2 <- dtmp %>% filter(GFP_RFP_b5 <= q2 & GFP_RFP_b5 > q1     & punc.RFP)
  data.q3 <- dtmp %>% filter(GFP_RFP_b5 <= q3 & GFP_RFP_b5 > q2     & punc.RFP)
  data.q4 <- dtmp %>% filter(GFP_RFP_b5 <= q4 & GFP_RFP_b5 > q3     & punc.RFP)
  data.q5 <- dtmp %>% filter(GFP_RFP_b5  > q4                        & punc.RFP)
  data <- list(data.q1, data.q2, data.q3, data.q4, data.q5)

  frac.punc.q1 <- nrow(data.q1) / sum(dtmp$GFP_RFP_b5 <  q1)
  frac.punc.q2 <- nrow(data.q2) / sum(dtmp$GFP_RFP_b5 <  q2 & dtmp$GFP_RFP_b5 > q1)
  frac.punc.q3 <- nrow(data.q3) / sum(dtmp$GFP_RFP_b5 <  q3 & dtmp$GFP_RFP_b5 > q2)
  frac.punc.q4 <- nrow(data.q4) / sum(dtmp$GFP_RFP_b5 <  q4 & dtmp$GFP_RFP_b5 > q3)
  frac.punc.q5 <- nrow(data.q5) / sum(dtmp$GFP_RFP_b5 >  q4)
  frac.punc.list <- list(round(frac.punc.q1,3), round(frac.punc.q2,3),
                         round(frac.punc.q3,3), round(frac.punc.q4,3),
                         round(frac.punc.q5,3))

  n.q <- sum(c(nrow(data.q1), nrow(data.q2), nrow(data.q3),
               nrow(data.q4), nrow(data.q5)) > 20)

  if (n.q > 0) {
    SUM       <- fig.data(n = n.q, FP = "RFP")
    plot.list[[I]] <- val.plot.full(n = n.q, FP = "RFP")
  } else {
    plot.list[[I]] <- val.plot.full(n = 0,    FP = "RFP")
  }
}

pdf(file.path(out.dir, "fig3C_TAF15_valency_scatter.pdf"), width = 6, height = 7)
print(plot.list[[1]])
dev.off()

cat("Wrote:", file.path(out.dir, "fig3C_TAF15_valency_scatter.pdf"), "\n")
