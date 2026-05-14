## debug_fig4D.R
## Diagnostic: trace fig4D's lookup-and-filter pipeline to find why
## D_score stays all NA. Run after sourcing fig4D once (or just standalone).

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("Fig4", "CD"))
data.KPlib <- dat$data
platedef   <- dat$platedef

cat("=== platedef ===\n")
cat("nrow:", nrow(platedef), "\n")
cat("col names:", paste(colnames(platedef), collapse=", "), "\n")
cat("unique GFP values (head 20):\n"); print(head(unique(platedef$GFP), 20))
cat("unique RFP values (head 20):\n"); print(head(unique(platedef$RFP), 20))

cat("\n=== data.KPlib ===\n")
cat("class:", class(data.KPlib), " length:", length(data.KPlib), "\n")
cat("first 10 names:\n"); print(head(names(data.KPlib), 10))
cat("nrow of first entry:", nrow(data.KPlib[[1]]), "\n")
cat("col names (first 20):", paste(head(colnames(data.KPlib[[1]]), 20), collapse=", "), "\n")

cat("\n=== filter test ===\n")
GFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43","crowder")
RFP.labels <- c("FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43",
                "Ddx4","DYRK3","Era","UBQ2","crowder")

cat("GFP.labels in platedef$GFP:\n")
print(GFP.labels %in% platedef$GFP)
cat("RFP.labels in platedef$RFP:\n")
print(RFP.labels %in% platedef$RFP)

valid.wells.no.data.filter <- platedef$WELL[
  platedef$GFP %in% GFP.labels &
  platedef$RFP %in% RFP.labels
]
cat("wells matching GFP+RFP filter:", length(valid.wells.no.data.filter), "\n")
cat("first 10:", paste(head(valid.wells.no.data.filter, 10), collapse=", "), "\n")

valid.wells <- intersect(valid.wells.no.data.filter, names(data.KPlib))
cat("wells also present in data.KPlib:", length(valid.wells), "\n")
cat("first 10:", paste(head(valid.wells, 10), collapse=", "), "\n")

if (length(valid.wells) == 0) {
  cat("\n!!! NO WELLS PASS BOTH FILTERS — STOPPING !!!\n")
  cat("platedef$WELL sample:\n"); print(head(platedef$WELL, 10))
  cat("names(data.KPlib) sample:\n"); print(head(names(data.KPlib), 10))
  return(invisible())
}

cat("\n=== walk one well end-to-end (the first valid one) ===\n")
W <- valid.wells[1]
G <- platedef$GFP[platedef$WELL == W][1]
R <- platedef$RFP[platedef$WELL == W][1]
cat("well:", W, " GFP:", G, " RFP:", R, "\n")

d.f <- data.KPlib[[W]]
cat("nrow before filter:", nrow(d.f), "\n")

d.f.filt <- d.f %>% dplyr::filter(inGFPnfoci >= 1 & inRFPnfoci >= 1 & f1_inGFP_toGFPmean > 150)
cat("nrow after filter (inGFPnfoci>=1 & inRFPnfoci>=1 & f1_inGFP_toGFPmean>150):",
    nrow(d.f.filt), "\n")

# Dissect the filter
cat("\nbreakdown:\n")
cat("  total rows:                              ", nrow(d.f), "\n")
cat("  inGFPnfoci >= 1:                          ", sum(as.numeric(d.f$inGFPnfoci) >= 1, na.rm = TRUE), "\n")
cat("  inRFPnfoci >= 1:                          ", sum(as.numeric(d.f$inRFPnfoci) >= 1, na.rm = TRUE), "\n")
cat("  both >= 1:                                ",
    sum(as.numeric(d.f$inGFPnfoci) >= 1 & as.numeric(d.f$inRFPnfoci) >= 1, na.rm = TRUE), "\n")
cat("  f1_inGFP_toGFPmean > 150:                 ",
    sum(as.numeric(d.f$f1_inGFP_toGFPmean) > 150, na.rm = TRUE), "\n")
cat("  range of f1_inGFP_toGFPmean:              ",
    paste(round(range(as.numeric(d.f$f1_inGFP_toGFPmean), na.rm = TRUE), 2), collapse=" - "), "\n")

if (nrow(d.f.filt) > 0) {
  cat("\nsample of f1 GFP foci coords (first 3):\n")
  print(head(d.f.filt[, c("f1_inGFPx", "f1_inGFPy", "f1_inRFPx", "f1_inRFPy")], 3))
}

cat("\n=== done ===\n")
