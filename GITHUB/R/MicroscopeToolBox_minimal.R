## MicroscopeToolBox_minimal.R
## Stripped-down version of the upstream MicroscopeToolBox.R from
## /Users/atarg/Documents/Levy_lab/scripts/RToolBox2/MicroscopeToolBox.R
## Only exports the helpers that the runners in scripts/ actually use.
##
## - filter.saturated() and correct.function() are NOT defined here because
##   identical copies live in functions_repository.R (the _w_sampling variant).
## - get.sNUM() is included for ad-hoc inspection of which microscope fields
##   (s_###) belong to a given well; not strictly required by any runner but
##   handy when debugging a platedef-vs-data mismatch.

get.sNUM <- function(design, plate = 1, well) {
  positions   <- design$WELLS[[plate]] %in% c(well)
  s.positions <- design$S.positions[[plate]][positions]
  return(as.numeric(s.positions))
}
