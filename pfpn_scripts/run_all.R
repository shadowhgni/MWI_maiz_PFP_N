# =============================================================================
# run_all.R: the whole PFP-N pipeline, in order
# -----------------------------------------------------------------------------
# Please, keep all scripts together in one folder 
# and run this script (run_all.R) with that folder as the working directory.
#
# Every output folder (pfpn_lookup/, pfpn_config/, pfpn_compiled/, pfpn_analysis/,
# pfpn_analysis_2/, covariates/) is created under `pfpn_root` if missing, so the
# pipeline runs from empty folders. Needs internet on the first run (World Bank
# DDI codebooks, terminag, Carob, GADM, covariates) and the LSMS zips whose paths
# are set in pfpn_compile.R.
#
# Before deleting pfpn_config/: it holds the curated tables  (unit weights,
# fertilizer N contents, crop names). They are rebuilt from the data and the
# corrections built into the scripts, but values you typed in by hand are lost.
# =============================================================================

# clean global environment
rm(list = ls())

# load core packages
library(tidyverse)
# Other packages used via pkg::fun (none)

pfpn_root <- ".."               # relative path to the scripts folder

options(
  pfpn.root    = pfpn_root,
  pfpn.refresh = "if_updated"      # "if_updated" (re-download only newer files), "never", "always"
)

steps <- c(
  "01.pfpn_variable_lookup.R",     # variable spec checked against each DDI -> pfpn_config/source_map.csv
  "02.pfpn_compile.R",             # LSMS + LCAS + Carob -> pfpn_compiled/pfpn_unified.rds (+ dictionary)
  "03.pfpn_analysis.R",            # data robustness, GAMs, variance components
  "04.pfpn_analysis_2.R"           # response curves, random forests + SHAP, report (Word + HTML)
)

missing <- steps[!file.exists(steps)]
if (length(missing) > 0) stop("Not found in ", getwd(), ": ", paste(missing, collapse = ", "))

# The loop runs inside a function, and each script in its own environment, so
# an rm(list = ls()) at the top of any script cannot remove what this loop needs.
run_steps <- function(steps) {
  for (step in steps) {
    t0 <- Sys.time()
    message("\n==================== ", step, " ====================")
    source(step, local = new.env(parent = globalenv()), echo = FALSE)
    message("---- ", step, " finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  }
  message("\nAll steps done. Report: ",
          file.path(getOption("pfpn.root", "."), "pfpn_analysis_2", "pfpn_report.docx"))
}
run_steps(steps)

# ==============================================================================
# END
# ==============================================================================