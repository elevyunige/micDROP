## fig4F_TDP43_specificity_heatmap.R
## Reproduces Fig 4F: TDP43 WT and three deletion variants (rows) x 11
## RFP IDRs (cols) colocalization heatmap.
## Source: pheatmap.R, run on the Fig4/F plate.
##
## Construct names match Fig4/F/platedef_idr.csv exactly (updated platedef
## uses "TDP43_331_343" in the GFP column, matching the published figure
## label "TDP43_Δ331-343").
##
## See fig4D_colocalization_heatmap.R for:
##   - the rationale behind the na.rm = TRUE inside both the rowwise
##     min() and the median();
##   - the convention that the "Valency" platedef label is the empty
##     2CLB-Dps control, displayed in the heatmap as "control";
##   - the rationale for using a permissive MIN_CELLS = 30 threshold
##     instead of the original 60;
##   - the removal of the brightness threshold (>150) from the filter:
##     the filter is just `inGFPnfoci >= 1 & inRFPnfoci >= 1`.

MIN_CELLS <- 30

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("Fig4", "F"))
data.KPlib <- dat$data
platedef   <- dat$platedef

# Note: this plate uses ER_a (with underscore), whereas Fig4/CD used Era.
GFP.labels <- c("TDP43_WT","TDP43_320_330","TDP43_331_343","TDP43_320_343")
RFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43",
                "Ddx4","DYRK3","ER_a","UBQ2","Valency")

GFP.display <- GFP.labels
RFP.display <- RFP.labels; RFP.display[RFP.display == "Valency"] <- "control"

D_score <- matrix(NA, nrow = length(GFP.labels), ncol = length(RFP.labels),
                  dimnames = list(GFP.labels, RFP.labels))

valid.wells <- platedef$WELL[
  platedef$GFP %in% GFP.labels &
  platedef$RFP %in% RFP.labels &
  platedef$WELL %in% names(data.KPlib)
]

n.ok <- 0L; n.too.small <- 0L; n.empty <- 0L
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

expected.wells <- platedef$WELL[
  platedef$GFP %in% GFP.labels & platedef$RFP %in% RFP.labels
]
missing.wells <- setdiff(expected.wells, names(data.KPlib))

cat(sprintf("[fig4F] %d ok, %d below MIN_CELLS=%d, %d empty, %d expected wells absent from data.KPlib\n",
            n.ok, n.too.small, MIN_CELLS, n.empty, length(missing.wells)))
if (length(skipped.too.small) > 0) {
  cat("[fig4F] wells nulled by MIN_CELLS threshold:\n")
  for (s in skipped.too.small) cat("  -", s, "\n")
}
if (length(missing.wells) > 0) {
  cat("[fig4F] wells absent from data.KPlib (matrix cells will be NA):\n")
  for (W in missing.wells) {
    G <- platedef$GFP[platedef$WELL == W][1]
    R <- platedef$RFP[platedef$WELL == W][1]
    cat(sprintf("  - %s  G=%-15s  R=%s\n", W, G, R))
  }
}

mode(D_score) <- "numeric"

rownames(D_score) <- GFP.display
colnames(D_score) <- RFP.display

pdf(file.path(out.dir, "fig4F_TDP43_specificity_heatmap.pdf"), width = 16, height = 6)
pheatmap(D_score,
         cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "black", na_col = "gray97",
         breaks = c(0, 0.4, 0.6, 1.5),
         color  = c("lightgoldenrod1", "orange1", "firebrick1"),
         main   = "TDP43 deletion variants - Foci Co-localization (median distance, um)",
         display_numbers = round(D_score, 2),
         fontsize = 9)
dev.off()

cat("Wrote:", file.path(out.dir, "fig4F_TDP43_specificity_heatmap.pdf"), "\n")
