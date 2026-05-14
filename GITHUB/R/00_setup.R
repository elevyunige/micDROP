## 00_setup.R
## Common setup sourced.
## Loads libraries, defines paths, and sources the plotting functions.
##
## EXPECTED USAGE (from a runner):
##   here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
##   source(file.path(here.repo.root, "R", "00_setup.R"))
##
## EDIT data.root BELOW TO POINT AT YOUR LOCAL DATA DIRECTORY.
## Expected layout under data.root/:
##   Fig1/         platedef_idr.csv  results/raw_list.RDS
##   Fig2/2D/      platedef_idr.csv  results/raw_list.RDS
##   Fig2/2E/      ...
##   Fig3/         ...
##   Fig4/CD/      ...
##   Fig4/F/       ...
##   FigS2/  FigS4/  FigS5/  FigS6/  FigS8/  FigS9/  FigS10/
##   FigS11/A/  FigS11/B/
##   CorrelationPlots/corr.RDS

# ---- USER CONFIG -----------------------------------------------------
data.root <- "/Users/..PATH../DATA/"
# ----------------------------------------------------------------------

if (!exists("here.repo.root")) {
  stop("here.repo.root is not set. Each runner must set it before sourcing 00_setup.R.")
}

out.dir <- file.path(here.repo.root, "output")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

options(scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(tidyverse)
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
  library(grid)
  library(reshape2)
  library(RColorBrewer)
})

source(file.path(here.repo.root, "R", "functions_repository.R"))
source(file.path(here.repo.root, "R", "MicroscopeToolBox_minimal.R"))

# Helper to load + correct a per-figure dataset.
load.figure.data <- function(subdir) {
  fig.dir <- file.path(data.root, subdir)
  if (!dir.exists(fig.dir)) {
    stop("Could not find data subdirectory: ", fig.dir,
         "\nEdit data.root in R/00_setup.R.")
  }
  rds  <- file.path(fig.dir, "results", "raw_list.RDS")
  pdef <- file.path(fig.dir, "platedef_idr.csv")
  if (!file.exists(rds))  stop("Missing raw_list.RDS: ", rds)
  if (!file.exists(pdef)) stop("Missing platedef_idr.csv: ", pdef)

  # raw_list.RDS comes in two shapes across the data folders shipped here:
  #   (a) FLAT - readRDS() returns a named list of data.frames keyed by
  #       well IDs (e.g. "A01", "E01").  This is the shape used by
  #       Fig2/2D, Fig3, ..., i.e. all packaged figure data.
  #   (b) WRAPPED - readRDS() returns a list of length 1 whose [[1]] is
  #       the named list of data.frames.  This was the shape assumed by
  #       the original R_scriptsrRedo/ scripts.
  # We accept either and unwrap to the flat list-of-data.frames.
  raw <- readRDS(rds)
  if (is.list(raw) && length(raw) > 0 && is.data.frame(raw[[1]])) {
    data.list <- raw
  } else if (is.list(raw) && length(raw) == 1 &&
             is.list(raw[[1]]) && length(raw[[1]]) > 0 &&
             is.data.frame(raw[[1]][[1]])) {
    data.list <- raw[[1]]
  } else {
    stop("Unexpected raw_list.RDS structure in ", rds,
         "\n  top-level class: ", paste(class(raw), collapse = ","),
         "\n  length: ", length(raw),
         "\n  class of first element: ",
         if (length(raw) >= 1) paste(class(raw[[1]]), collapse = ",") else "NA")
  }

  # The verbatim helpers in R/functions_repository.R iterate by integer
  # index and assume each list element is a data.frame.  Drop NULL or
  # non-data.frame slots (e.g. BLANK wells stored as NULL) before they
  # hit the helpers.  Names of kept slots are preserved, so later
  # look-ups of the form data.list[[W]] still work.
  is.df <- vapply(data.list,
                  function(x) is.data.frame(x) || is.matrix(x),
                  logical(1))
  data.list <- data.list[is.df]

  data.list <- filter.saturated(data.list, "GFP_int_b0")
  data.list <- filter.saturated(data.list, "RFP_int_b0")
  data.list <- correct.function(data.list, "GFP_int_b5",         "RFP_int_b5")
  data.list <- correct.function(data.list, "GFP_int_b0",         "RFP_int_b0")
  data.list <- correct.function(data.list, "GFP_int_mean",       "RFP_int_mean")
  data.list <- correct.function(data.list, "f1_inGFP_toGFPmean", "f1_inGFP_toRFPmean")

  pdef.df <- read.csv(pdef)
  # Standardize column names: ensure WELL, BF, GFP, RFP are upper-case.
  required <- c("well", "bf", "gfp", "rfp")
  hits <- match(required, tolower(colnames(pdef.df)))
  if (anyNA(hits)) {
    stop("Platedef ", pdef, " is missing one of: WELL, BF, GFP, RFP.")
  }
  colnames(pdef.df)[hits] <- c("WELL", "BF", "GFP", "RFP")

  list(data = data.list, platedef = pdef.df)
}
