## fig4D_colocalization_heatmap.R
## Reproduces Fig 4D: 7 GFP-IDRs x 11 RFP-IDRs colocalization heatmap.
## Source: pheatmap.R, run on the Fig4/CD plate.
##
## Implementation notes:
##
## 1. Brightness filter.
##    Cells are kept only if `inGFPnfoci >= 1 & inRFPnfoci >= 1 &
##    f1_inGFP_toGFPmean > 150`. The 150 threshold rejects cells whose
##    GFP focus brightness is below the imaging noise floor — without
##    it, false-positive single-pixel detections contaminate the
##    median distance. This filter is what produces the published
##    Fig 4D, so we keep it here.
##
## 2. Minimum cells per well.
##    MIN_CELLS = 60, matching the original pheatmap.R threshold. A
##    well with fewer than 60 cells passing the brightness filter is
##    reported as NA. The threshold is checked twice: once on the
##    per-well row count after filtering, and once on the count of
##    cells that yield a finite distance after NA propagation.
##
## 3. Distance handling for cells with < 2 foci.
##    For every cell we compute four candidate distances (between the
##    two possible GFP foci and the two possible RFP foci).  When a
##    channel has only one focus, the corresponding f2_* columns are
##    NA. `valid_pair_distance` returns NA whenever any of the four
##    coordinates is non-finite, and we additionally invalidate
##    distance2 / distance3 / distance4 explicitly when they involve
##    a non-existent second focus (by inspecting inGFPnfoci /
##    inRFPnfoci).  `row_min_na` then takes the minimum of the four
##    candidates while gracefully handling rows where all four are NA.
##
## 4. Control row & column.
##    The "Valency" platedef label corresponds to the empty 2CLB-Dps
##    scaffold fused to a fluorophore but with no IDR. It is the
##    "no IDR" control and is displayed as "control" in the heatmap.
##
## 5. Single-page output (raw distance only).
##    The original pheatmap.R produced 3 pages: raw distance,
##    GFP-normalized (each row / control row), RFP-normalized (each
##    col / control col). Two issues with that:
##      a) pheatmap() does not call grid.newpage() automatically in
##         a pdf() device, so the three calls overlay onto a single
##         page and only the last (RFP-normalized) is visible.
##      b) The normalized pages divide by the control row/col, which
##         propagates any single NA in the control into an entire
##         row/col of NAs — producing the "lots of NA" artifact even
##         when the underlying raw matrix is well populated.
##    The published Fig 4D shows raw distances anyway, so this runner
##    only emits the raw heatmap. The raw matrix is also written
##    alongside as a CSV for inspection.

MIN_CELLS  <- 60
PIXEL_SIZE <- 0.108

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

dat        <- load.figure.data(file.path("Fig4", "CD"))
data.KPlib <- dat$data
platedef   <- dat$platedef

GFP.list.full <- c("FUS", "hnRNPA1", "HspB8", "RBM14", "TAF15", "TDP43", "Valency")
RFP.list.full <- c("FUS", "hnRNPA1", "HspB8", "RBM14", "TAF15", "TDP43",
                   "Ddx4", "DYRK3", "Era", "UBQ2", "Valency")

GFP.display <- GFP.list.full; GFP.display[GFP.display == "Valency"] <- "control"
RFP.display <- RFP.list.full; RFP.display[RFP.display == "Valency"] <- "control"

# ---- distance helpers ----------------------------------------------

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

distance_xy <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

valid_pair_distance <- function(gx, gy, rx, ry) {
  d  <- distance_xy(gx, gy, rx, ry)
  ok <- is.finite(gx) & is.finite(gy) & is.finite(rx) & is.finite(ry)
  d[!ok] <- NA_real_
  d
}

row_min_na <- function(mat) {
  out <- apply(mat, 1, function(z) {
    if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)
  })
  as.numeric(out)
}

# Returns a list with $median (the well's D_score) and $n_pass (the
# number of cells that yielded a finite per-cell distance), so the
# caller can report exactly why a well was nulled.
calc_D_score_for_well <- function(df, pixel_size = PIXEL_SIZE) {
  if (is.null(df) || nrow(df) == 0) {
    return(list(median = NA_real_, n_pass = 0L, n_after_filter = 0L))
  }

  needed <- c("inGFPnfoci", "inRFPnfoci", "f1_inGFP_toGFPmean",
              "f1_inGFPx", "f1_inGFPy", "f1_inRFPx", "f1_inRFPy",
              "f2_inGFPx", "f2_inGFPy", "f2_inRFPx", "f2_inRFPy")
  if (!all(needed %in% colnames(df))) {
    return(list(median = NA_real_, n_pass = 0L, n_after_filter = 0L))
  }

  df <- df %>%
    mutate(across(all_of(needed), safe_numeric)) %>%
    filter(inGFPnfoci >= 1,
           inRFPnfoci >= 1,
           f1_inGFP_toGFPmean > 150)

  n_after_filter <- nrow(df)
  if (n_after_filter <= MIN_CELLS) {
    return(list(median = NA_real_, n_pass = 0L,
                n_after_filter = n_after_filter))
  }

  df$distance1 <- valid_pair_distance(df$f1_inGFPx, df$f1_inGFPy,
                                      df$f1_inRFPx, df$f1_inRFPy)
  df$distance2 <- valid_pair_distance(df$f2_inGFPx, df$f2_inGFPy,
                                      df$f2_inRFPx, df$f2_inRFPy)
  df$distance3 <- valid_pair_distance(df$f1_inGFPx, df$f1_inGFPy,
                                      df$f2_inRFPx, df$f2_inRFPy)
  df$distance4 <- valid_pair_distance(df$f2_inGFPx, df$f2_inGFPy,
                                      df$f1_inRFPx, df$f1_inRFPy)

  # Explicitly invalidate distances involving nonexistent second foci.
  df$distance2[df$inGFPnfoci < 2 | df$inRFPnfoci < 2] <- NA_real_
  df$distance3[df$inRFPnfoci < 2] <- NA_real_
  df$distance4[df$inGFPnfoci < 2] <- NA_real_

  df$distance <- row_min_na(df[, c("distance1", "distance2", "distance3", "distance4")])
  df$distance <- df$distance * pixel_size

  n_pass <- sum(is.finite(df$distance))
  if (n_pass <= MIN_CELLS) {
    return(list(median = NA_real_, n_pass = n_pass,
                n_after_filter = n_after_filter))
  }

  list(median = median(df$distance, na.rm = TRUE),
       n_pass = n_pass,
       n_after_filter = n_after_filter)
}

# ---- fill the matrix -----------------------------------------------

D_score <- matrix(NA_real_, nrow = length(GFP.list.full), ncol = length(RFP.list.full),
                  dimnames = list(GFP.list.full, RFP.list.full))

# Per-well diagnostics matrix (cells passing the brightness filter,
# cells passing all the way through to a finite distance).
N_after_filter <- matrix(NA_integer_, nrow = nrow(D_score), ncol = ncol(D_score),
                         dimnames = dimnames(D_score))
N_pass         <- matrix(NA_integer_, nrow = nrow(D_score), ncol = ncol(D_score),
                         dimnames = dimnames(D_score))

n.ok <- 0L; n.empty <- 0L; n.absent <- 0L; n.below <- 0L
absent.cells <- character(0)
below.cells  <- character(0)

for (I in seq_along(GFP.list.full)) {
  for (J in seq_along(RFP.list.full)) {
    G <- GFP.list.full[I]
    R <- RFP.list.full[J]
    W <- platedef$WELL[which(platedef$GFP == G & platedef$RFP == R)]

    if (length(W) != 1 || !(W %in% names(data.KPlib))) {
      D_score[I, J] <- NA_real_
      n.absent <- n.absent + 1L
      absent.cells <- c(absent.cells, sprintf("G=%s R=%s", G, R))
      next
    }

    DAT <- data.KPlib[[W]]
    if (is.null(DAT) || nrow(DAT) == 0) {
      D_score[I, J] <- NA_real_
      n.empty <- n.empty + 1L
      next
    }

    res <- calc_D_score_for_well(DAT, pixel_size = PIXEL_SIZE)
    D_score[I, J]        <- res$median
    N_after_filter[I, J] <- as.integer(res$n_after_filter)
    N_pass[I, J]         <- as.integer(res$n_pass)
    if (is.finite(res$median)) {
      n.ok <- n.ok + 1L
    } else {
      n.below <- n.below + 1L
      below.cells <- c(below.cells,
                       sprintf("%s (G=%s R=%s, n_filt=%d, n_pass=%d)",
                               W, G, R, res$n_after_filter, res$n_pass))
    }
  }
}

mode(D_score) <- "numeric"

cat(sprintf("[fig4D] %d ok, %d below MIN_CELLS=%d (after >150 brightness filter), %d empty wells, %d expected cells absent from data.KPlib\n",
            n.ok, n.below, MIN_CELLS, n.empty, n.absent))
if (length(below.cells) > 0) {
  cat("[fig4D] wells nulled by MIN_CELLS / brightness threshold:\n")
  for (s in below.cells) cat("  -", s, "\n")
}
if (length(absent.cells) > 0) {
  cat("[fig4D] absent / un-mappable matrix cells:\n")
  for (s in absent.cells) cat("  -", s, "\n")
}

# ---- write the raw matrix as CSV (for inspection) ------------------

D_for_csv <- D_score
rownames(D_for_csv) <- GFP.display
colnames(D_for_csv) <- RFP.display
write.csv(round(D_for_csv, 3),
          file.path(out.dir, "fig4D_colocalization_matrix.csv"),
          row.names = TRUE, na = "NA")

# Diagnostic per-well counts (cells passing the >150 filter, cells
# yielding a finite distance), to help judge whether MIN_CELLS=60 is
# appropriate for this dataset.
N_for_csv <- N_after_filter
rownames(N_for_csv) <- GFP.display; colnames(N_for_csv) <- RFP.display
write.csv(N_for_csv,
          file.path(out.dir, "fig4D_colocalization_n_after_filter.csv"),
          row.names = TRUE, na = "NA")

P_for_csv <- N_pass
rownames(P_for_csv) <- GFP.display; colnames(P_for_csv) <- RFP.display
write.csv(P_for_csv,
          file.path(out.dir, "fig4D_colocalization_n_pass.csv"),
          row.names = TRUE, na = "NA")

# ---- plot (single page, raw distance only) -------------------------

M <- D_score
rownames(M) <- GFP.display
colnames(M) <- RFP.display

pdf(file.path(out.dir, "fig4D_colocalization_heatmap.pdf"), width = 16, height = 6)
pheatmap(M,
         cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "black", na_col = "gray97",
         breaks = c(0, 0.4, 0.6, 1.5),
         color  = c("lightgoldenrod1", "orange1", "firebrick1"),
         main   = "Foci Co-localization\nMedian Distance (um)",
         display_numbers = round(M, 2),
         fontsize = 9)
dev.off()

cat("Wrote:", file.path(out.dir, "fig4D_colocalization_heatmap.pdf"), "\n")
cat("Wrote:", file.path(out.dir, "fig4D_colocalization_matrix.csv"),  "\n")
cat("Wrote:", file.path(out.dir, "fig4D_colocalization_n_after_filter.csv"), "\n")
cat("Wrote:", file.path(out.dir, "fig4D_colocalization_n_pass.csv"),  "\n")
