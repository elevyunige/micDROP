## 01_run.R
## Sources every working runner in scripts/ in figure order, captures
## errors per runner so one failure doesn't halt the rest, and prints a
## summary table at the end.
##
## Usage from an R session at the repository root:
##   source("R/01_run.R")
##
## All output PDFs are written to output/, exactly as if each runner had
## been sourced individually.
##
## NOTE on path resolution.  Each runner's first non-comment line is
##   here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
## which, when the runner is sourced from this dispatcher, resolves
## sys.frame(1)$ofile to *this* file (R/01_run.R) rather than the runner
## itself.  That's intentionally fine: dirname("R/01_run.R")/.. and
## dirname("scripts/figXX.R")/.. both point at the repo root because
## R/ and scripts/ are sibling directories.  If you ever move this file
## out of R/, that coincidence breaks — keep 01_run.R alongside 00_setup.R.
##
## NOTE on R scoping.  Runners are sourced into the global environment
## (no `local = new.env()`), because the helpers in functions_repository.R
## (val.plot.full, fig.data, medmax.plot, ...) look up free variables
## (`data`, `frac.punc.q1`, `dtmp`, etc.) via lexical scoping in their
## enclosing environment.  Sourcing 00_setup.R into globalenv puts those
## helpers' enclosure at globalenv, so the runner's variables must also
## be in globalenv for the helpers to find them.  Variable bleed between
## runners is harmless because every runner overwrites the variables it
## uses at the top of its body.
##
## NOTE on loop variable hygiene.  Because runners share globalenv with
## this dispatcher, any name a runner uses for its own for-loop index
## clobbers the dispatcher's index too.  We therefore use distinctive
## names (`run.i`, `fail.i`) for the dispatcher's loops; the runners
## should avoid these specific names.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
scripts.dir <- file.path(here.repo.root, "scripts")
if (!dir.exists(scripts.dir)) {
  stop("Could not find scripts/ at ", scripts.dir)
}

# ---- Runner discovery -------------------------------------------------
all.runners <- list.files(scripts.dir, pattern = "\\.R$", full.names = TRUE)

# Drop ad-hoc debug/diagnostic scripts (debug_*.R, scratch_*.R, etc.) so
# the dispatcher only ever runs proper figure runners. To run a debug
# script, source it directly: source("scripts/debug_fig4D.R").
all.runners <- all.runners[
  !grepl("^(debug|scratch|tmp|wip)_", basename(all.runners), ignore.case = TRUE)
]

# Skip deprecated stubs (their first non-comment line is `stop(...)`).
is.stub <- function(path) {
  txt <- readLines(path, warn = FALSE)
  txt <- txt[!grepl("^\\s*(#|$)", txt)]
  length(txt) > 0 && grepl("^\\s*stop\\(", txt[1])
}
all.runners <- all.runners[!vapply(all.runners, is.stub, logical(1))]

# Sort into canonical figure order: main figs (fig2…fig4) first by
# numeric+suffix, then supplementary (figS…) by numeric+suffix. The
# fallback uses list() — not c() — because c(2L, 999L, "z") would
# coerce all three to character, breaking vapply(..., integer(1)) below.
sort.key <- function(path) {
  fname <- basename(path)
  m <- regmatches(fname, regexec("^fig(S?)(\\d+)([A-Za-z]*)_", fname))[[1]]
  if (length(m) < 4) return(list(2L, 999L, "z"))
  main.or.supp <- if (m[2] == "") 0L else 1L
  fig.num      <- as.integer(m[3])
  fig.suffix   <- m[4]
  list(main.or.supp, fig.num, fig.suffix)
}
keys <- lapply(all.runners, sort.key)
ord <- order(vapply(keys, `[[`, integer(1), 1),
             vapply(keys, `[[`, integer(1), 2),
             vapply(keys, `[[`, character(1), 3))
runners <- all.runners[ord]

cat("Found", length(runners), "runner(s) to execute, in this order:\n")
cat(paste0("  - ", basename(runners), collapse = "\n"), "\n\n")

# ---- Dispatch loop ----------------------------------------------------
results <- data.frame(
  script  = basename(runners),
  status  = NA_character_,
  seconds = NA_real_,
  message = NA_character_,
  stringsAsFactors = FALSE
)

for (run.i in seq_along(runners)) {
  cat(strrep("=", 70), "\n", sep = "")
  cat("Running ", basename(runners[run.i]),
      "  (", run.i, "/", length(runners), ")\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")

  t0  <- Sys.time()
  err <- NULL
  ok <- tryCatch({
    # NOTE: source() into the global environment (no `local =`) — see
    # the scoping note at the top of this file.
    source(runners[run.i])
    TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
    FALSE
  })

  # Re-fetch the dispatcher's loop index after source(): a runner may
  # have shadowed it via its own for-loop. We recover by indexing on
  # status==NA (the as-yet-unwritten row).
  cur.row <- which(is.na(results$status))[1]
  if (is.na(cur.row)) cur.row <- length(runners)  # safety fallback

  results$seconds[cur.row] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results$status[cur.row]  <- if (ok) "OK" else "FAILED"
  results$message[cur.row] <- if (ok) "" else err

  cat("\n>>> ", results$status[cur.row], " in ",
      results$seconds[cur.row], "s\n\n", sep = "")
}

# ---- Summary ----------------------------------------------------------
cat(strrep("=", 70), "\n", sep = "")
cat("Run summary\n")
cat(strrep("=", 70), "\n", sep = "")
print(results, row.names = FALSE)

n.ok   <- sum(results$status == "OK",     na.rm = TRUE)
n.fail <- sum(results$status == "FAILED", na.rm = TRUE)
n.miss <- sum(is.na(results$status))
cat(sprintf("\n%d/%d runners completed successfully (%d failed, %d not recorded).\n",
            n.ok, length(runners), n.fail, n.miss))

if (n.fail > 0) {
  cat("\nFailures:\n")
  for (fail.i in which(results$status == "FAILED")) {
    cat("  - ", results$script[fail.i], ": ", results$message[fail.i], "\n", sep = "")
  }
}

invisible(results)
