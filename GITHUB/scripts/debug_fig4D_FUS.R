## debug_fig4D_FUS.R
## Print exhaustive details for every well in Fig4/CD where FUS is on
## either channel - i.e. GFP=FUS (the F row) or RFP=FUS (column 10) -
## paired with any of the IDRs in the heatmap.
##
## Filter is `inGFPnfoci >= 1 & inRFPnfoci >= 1` (no brightness threshold).
##
## Run as:
##   source("scripts/debug_fig4D_FUS.R")

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("Fig4", "CD"))
data.KPlib <- dat$data
platedef   <- dat$platedef

GFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43","Valency")
RFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43",
                "Ddx4","DYRK3","Era","UBQ2","Valency")

print.well.details <- function(W) {
  if (!(W %in% names(data.KPlib))) {
    cat(sprintf("%s   (well not in data.KPlib - SKIPPING)\n\n", W))
    return(invisible())
  }
  G <- platedef$GFP[platedef$WELL == W][1]
  R <- platedef$RFP[platedef$WELL == W][1]

  if (G == "FUS") {
    side <- "GFP=FUS"; partner <- R; partner.side <- "RFP"
  } else {
    side <- "RFP=FUS"; partner <- G; partner.side <- "GFP"
  }

  d <- data.KPlib[[W]]
  cat(strrep("=", 72), "\n", sep = "")
  cat(sprintf("%s   %s   partner: %s=%s\n", W, side, partner.side, partner))
  cat(strrep("-", 72), "\n", sep = "")

  raw.n <- nrow(d)

  for (col in c("inGFPnfoci","inRFPnfoci",
                "f1_inGFPx","f1_inGFPy","f1_inRFPx","f1_inRFPy",
                "f2_inGFPx","f2_inGFPy","f2_inRFPx","f2_inRFPy")) {
    if (col %in% names(d)) d[[col]] <- as.numeric(d[[col]])
  }

  n.gfp.foci  <- sum(d$inGFPnfoci >= 1, na.rm = TRUE)
  n.rfp.foci  <- sum(d$inRFPnfoci >= 1, na.rm = TRUE)
  n.passes    <- sum(d$inGFPnfoci >= 1 & d$inRFPnfoci >= 1, na.rm = TRUE)

  cat(sprintf("  raw cells:                                   %5d\n", raw.n))
  cat(sprintf("  cells with inGFPnfoci >= 1:                  %5d\n", n.gfp.foci))
  cat(sprintf("  cells with inRFPnfoci >= 1:                  %5d\n", n.rfp.foci))
  cat(sprintf("  cells passing the full filter (both >= 1):   %5d\n", n.passes))

  d.f <- d %>% dplyr::filter(inGFPnfoci >= 1 & inRFPnfoci >= 1)
  if (nrow(d.f) == 0) {
    cat("  (no cells survive the filter - cannot compute distances)\n\n")
    return(invisible())
  }

  one.foci <- d.f$inGFPnfoci == 1 & d.f$inRFPnfoci == 1
  cat(sprintf("  among passers, one-focus  (inG=inR=1):       %5d\n", sum(one.foci)))
  cat(sprintf("  among passers, multi-focus:                  %5d\n", sum(!one.foci)))

  d.f$f2_inRFPx[one.foci] <- 100000
  d.f$f2_inRFPy[one.foci] <- 100000
  d.f$distance1 <- sqrt((d.f$f1_inGFPx - d.f$f1_inRFPx)^2 + (d.f$f1_inGFPy - d.f$f1_inRFPy)^2)
  d.f$distance2 <- sqrt((d.f$f2_inGFPx - d.f$f2_inRFPx)^2 + (d.f$f2_inGFPy - d.f$f2_inRFPy)^2)
  d.f$distance3 <- sqrt((d.f$f1_inGFPx - d.f$f2_inRFPx)^2 + (d.f$f1_inGFPy - d.f$f2_inRFPy)^2)
  d.f$distance4 <- sqrt((d.f$f2_inGFPx - d.f$f1_inRFPx)^2 + (d.f$f2_inGFPy - d.f$f1_inRFPy)^2)
  d.f$distance  <- apply(d.f[, c("distance1","distance2","distance3","distance4")],
                         1, function(r) min(r, na.rm = TRUE)) * 0.108

  med <- median(d.f$distance, na.rm = TRUE)
  qd  <- quantile(d.f$distance,
                  probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  cat(sprintf("  median distance:                             %6.3f um\n", med))
  cat(sprintf("  distance q10/25/50/75/90: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              qd[1], qd[2], qd[3], qd[4], qd[5]))
  if (nrow(d.f) >= 30) {
    cat("  in heatmap (MIN_CELLS=30)?: YES\n")
  } else {
    cat(sprintf("  in heatmap (MIN_CELLS=30)?: NO  (only %d cells - cell will be NA)\n",
                nrow(d.f)))
  }
  cat("\n")
}

cat("################################################################\n")
cat("##                       FUS as GFP                            ##\n")
cat("################################################################\n\n")
gfp.fus.wells <- sort(platedef$WELL[platedef$GFP == "FUS" &
                                     platedef$RFP %in% RFP.labels])
for (W in gfp.fus.wells) print.well.details(W)

cat("################################################################\n")
cat("##                       FUS as RFP                            ##\n")
cat("################################################################\n\n")
rfp.fus.wells <- sort(platedef$WELL[platedef$RFP == "FUS" &
                                     platedef$GFP %in% GFP.labels])
for (W in rfp.fus.wells) print.well.details(W)

cat("################################################################\n")
cat("## Side-by-side: FUS-GFP vs FUS-RFP for each shared partner    ##\n")
cat("################################################################\n\n")

shared <- intersect(GFP.labels, RFP.labels)
shared <- setdiff(shared, c("FUS", "Valency"))

compute.median <- function(W) {
  if (!(W %in% names(data.KPlib))) return(NA_real_)
  d.f <- data.KPlib[[W]]
  if (nrow(d.f) == 0) return(NA_real_)
  for (col in c("f1_inGFPx","f1_inGFPy","f1_inRFPx","f1_inRFPy",
                "f2_inGFPx","f2_inGFPy","f2_inRFPx","f2_inRFPy",
                "inGFPnfoci","inRFPnfoci")) {
    if (col %in% names(d.f)) d.f[[col]] <- as.numeric(d.f[[col]])
  }
  d.f <- d.f %>% dplyr::filter(inGFPnfoci >= 1 & inRFPnfoci >= 1)
  if (nrow(d.f) < 30) return(NA_real_)
  one.foci <- d.f$inGFPnfoci == 1 & d.f$inRFPnfoci == 1
  d.f$f2_inRFPx[one.foci] <- 100000
  d.f$f2_inRFPy[one.foci] <- 100000
  d1 <- sqrt((d.f$f1_inGFPx - d.f$f1_inRFPx)^2 + (d.f$f1_inGFPy - d.f$f1_inRFPy)^2)
  d2 <- sqrt((d.f$f2_inGFPx - d.f$f2_inRFPx)^2 + (d.f$f2_inGFPy - d.f$f2_inRFPy)^2)
  d3 <- sqrt((d.f$f1_inGFPx - d.f$f2_inRFPx)^2 + (d.f$f1_inGFPy - d.f$f2_inRFPy)^2)
  d4 <- sqrt((d.f$f2_inGFPx - d.f$f1_inRFPx)^2 + (d.f$f2_inGFPy - d.f$f1_inRFPy)^2)
  dd <- pmin(d1, d2, d3, d4, na.rm = TRUE) * 0.108
  median(dd, na.rm = TRUE)
}

cat(sprintf("  %-10s  %-22s  %-22s\n",
            "partner", "FUS-GFP x partner-RFP", "partner-GFP x FUS-RFP"))
cat(sprintf("  %-10s  %-22s  %-22s\n",
            strrep("-", 10), strrep("-", 22), strrep("-", 22)))
for (P in shared) {
  W1 <- platedef$WELL[platedef$GFP == "FUS" & platedef$RFP == P][1]
  W2 <- platedef$WELL[platedef$GFP == P     & platedef$RFP == "FUS"][1]
  m1 <- if (!is.na(W1)) compute.median(W1) else NA_real_
  m2 <- if (!is.na(W2)) compute.median(W2) else NA_real_
  s1 <- if (is.na(m1)) "NA" else sprintf("%.3f um", m1)
  s2 <- if (is.na(m2)) "NA" else sprintf("%.3f um", m2)
  cat(sprintf("  %-10s  %-22s  %-22s\n", P, s1, s2))
}
W.self <- platedef$WELL[platedef$GFP == "FUS" & platedef$RFP == "FUS"][1]
m.self <- if (!is.na(W.self)) compute.median(W.self) else NA_real_
cat(sprintf("  %-10s  %-22s\n", "FUS (self)",
            if (is.na(m.self)) "NA" else sprintf("%.3f um", m.self)))
