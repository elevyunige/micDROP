## fig4D_colocalization_heatmap.R
## Reproduces Fig 4D: 7 GFP-IDRs x 11 RFP-IDRs colocalization heatmap.
## Source: pheatmap.R, run on the Fig4/CD plate.
##
## Implementation notes:
##
## 1. Colour scheme.
##    The script's hard-coded colour palette is c("lightgoldenrod1",
##    "orange1","firebrick1") — yellow -> orange -> red.  The published
##    Fig 4D uses a yellow -> pink -> magenta gradient.  Per author
##    decision (2026-04), this runner keeps the original yellow -> orange
##    -> red palette; the colour mismatch with the published figure is
##    intentional.
##
## 2. One-focus cells in the distance calculation.
##    For cells with exactly one focus on each channel, the f2 columns
##    are NA.  The original script masks ONLY the f2_inRFP columns to
##    100000 to push distance3 (= f1_GFP <-> f2_RFP) far away from the
##    minimum, but leaves f2_inGFP at NA.  Without na.rm on the rowwise
##    min and on the median, distance2 and distance4 propagate NA into
##    every one-focus row, and the median collapses to NA for every well.
##    Adding na.rm = TRUE to both calls preserves the intended semantics
##    (distance1 wins for one-focus cells; distance2..4 ignored when
##    they're undefined).
##
## 3. Control row & column.
##    The published Fig 4D's "control" row and column correspond to the
##    "Valency" platedef label — the empty 2CLB-Dps scaffold fused to
##    a fluorophore but without any IDR.  The "crowder" platedef label
##    (3IQ1-Dps based) is a different construct and not used in Fig 4D.
##
## 4. Filter.
##    Just `inGFPnfoci >= 1 & inRFPnfoci >= 1` — keep cells with at
##    least one focus on each channel, no brightness threshold.
##
## 5. Minimum cells per well.
##    MIN_CELLS = 30. The original pheatmap.R used 60 (line `if (nrow
##    (...) > 60)`); 30 is more permissive without admitting noisy
##    wells.
##
## 6. Normalized matrices dropped.
##    The original pheatmap.R produced 3 pages: raw distance, GFP-
##    normalized (each row divided by the control row), RFP-normalized
##    (each col divided by the control col).  Those normalizations
##    force the control row/col to be exactly 1.0 by construction
##    (control divided by itself), which is misleading.  The published
##    Fig 4D shows raw distances anyway, so this runner only emits the
##    raw heatmap.

stop("DEPRECATED: superseded by fig4D_colocalization_heatmap.R (kept for reference only)")

MIN_CELLS <- 30

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("Fig4", "CD"))
data.KPlib <- dat$data
platedef   <- dat$platedef

# Matrix axis labels in display order. The 7th/11th label is "Valency"
# in the platedef, displayed as "control" in the heatmap.
GFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43","Valency")
RFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43",
                "Ddx4","DYRK3","Era","UBQ2","Valency")

GFP.display <- GFP.labels; GFP.display[GFP.display == "Valency"] <- "control"
RFP.display <- RFP.labels; RFP.display[RFP.display == "Valency"] <- "control"

D_score <- matrix(NA, nrow = length(GFP.labels), ncol = length(RFP.labels),
                  dimnames = list(GFP.labels, RFP.labels))

valid.wells <- platedef$WELL[
  platedef$GFP %in% GFP.labels &
  platedef$RFP %in% RFP.labels &
  platedef$WELL %in% names(data.KPlib)
]

n.empty <- 0L; n.too.small <- 0L; n.ok <- 0L
skipped.too.small <- character(0)

for (W in valid.wells) {
  G <- platedef$GFP[platedef$WELL == W][1]
  R <- platedef$RFP[platedef$WELL == W][1]
  d.f <- data.KPlib[[W]]
  if (is.null(d.f) || nrow(d.f) == 0) { n.empty <- n.empty + 1L; next }

  d.f <- d.f %>% filter(inGFPnfoci >= 1 & inRFPnfoci >= 1)
  if (nrow(d.f) < MIN_CELLS) {
    n.too.small <- n.too.small + 1L
    skipped.too.small <- c(skipped.too.small,
                            sprintf("%s (G=%s R=%s, n=%d)", W, G, R, nrow(d.f)))
    next
  }

  for (col in c("f1_inGFPx","f1_inGFPy","f1_inRFPx","f1_inRFPy",
                "f2_inGFPx","f2_inGFPy","f2_inRFPx","f2_inRFPy")) {
    d.f[[col]] <- as.numeric(d.f[[col]])
  }
  one.foci <- d.f$inGFPnfoci == 1 & d.f$inRFPnfoci == 1
  d.f$f2_inRFPx[one.foci] <- 100000
  d.f$f2_inRFPy[one.foci] <- 100000

  d.f$distance1 <- sqrt((d.f$f1_inGFPx - d.f$f1_inRFPx)^2 + (d.f$f1_inGFPy - d.f$f1_inRFPy)^2)
  d.f$distance2 <- sqrt((d.f$f2_inGFPx - d.f$f2_inRFPx)^2 + (d.f$f2_inGFPy - d.f$f2_inRFPy)^2)
  d.f$distance3 <- sqrt((d.f$f1_inGFPx - d.f$f2_inRFPx)^2 + (d.f$f1_inGFPy - d.f$f2_inRFPy)^2)
  d.f$distance4 <- sqrt((d.f$f2_inGFPx - d.f$f1_inRFPx)^2 + (d.f$f2_inGFPy - d.f$f1_inRFPy)^2)
  d.f$distance  <- apply(d.f[, c("distance1","distance2","distance3","distance4")],
                         1, function(r) min(r, na.rm = TRUE)) * 0.108

  D_score[G, R] <- median(d.f$distance, na.rm = TRUE)
  n.ok <- n.ok + 1L
}

# Wells expected by the matrix but completely absent from data.KPlib
expected.wells <- platedef$WELL[
  platedef$GFP %in% GFP.labels & platedef$RFP %in% RFP.labels
]
missing.wells <- setdiff(expected.wells, names(data.KPlib))

cat(sprintf("[fig4D] %d ok, %d below MIN_CELLS=%d, %d empty, %d expected wells absent from data.KPlib\n",
            n.ok, n.too.small, MIN_CELLS, n.empty, length(missing.wells)))
if (length(skipped.too.small) > 0) {
  cat("[fig4D] wells nulled by MIN_CELLS threshold:\n")
  for (s in skipped.too.small) cat("  -", s, "\n")
}
if (length(missing.wells) > 0) {
  cat("[fig4D] wells absent from data.KPlib (matrix cells will be NA):\n")
  for (W in missing.wells) {
    G <- platedef$GFP[platedef$WELL == W][1]
    R <- platedef$RFP[platedef$WELL == W][1]
    cat(sprintf("  - %s  G=%-8s  R=%s\n", W, G, R))
  }
}

mode(D_score) <- "numeric"
rownames(D_score) <- GFP.display
colnames(D_score) <- RFP.display

pdf(file.path(out.dir, "fig4D_colocalization_heatmap.pdf"), width = 16, height = 6)
pheatmap(D_score,
         cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "black", na_col = "gray97",
         breaks = c(0, 0.4, 0.6, 1.5),
         color  = c("lightgoldenrod1", "orange1", "firebrick1"),
         main   = "Foci Co-localization\nMedian Distance (um)",
         display_numbers = round(D_score, 2),
         fontsize = 9)
dev.off()

cat("Wrote:", file.path(out.dir, "fig4D_colocalization_heatmap.pdf"), "\n")
