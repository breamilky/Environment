################################################################################
# 00_Bridge_RDS.R
#
# PURPOSE: Fix the filename mismatch that currently stops four scripts dead.
#
# THE PROBLEM
#   0202_first.R writes RDS objects with a "_2" suffix:
#       merged_daily_pm25_2.rds, first_stage_pm25_2.rds, results_pm25_2.rds ...
#   But 01_Feb2.R, 01_Feb2_2.R, 01_Feb2_3.R and 0202_second.R all read the
#   UNSUFFIXED names, and 01_Feb2_2.R / 01_Feb2_4_PM10.R additionally look for
#   "script02_all_results.RData".
#   Nothing in the 24-script codebase writes either the unsuffixed names or the
#   .RData file, so:
#       01_Feb2_2.R      -> stop("No saved data found. Run Script 02 first.")
#       01_Feb2_4_PM10.R -> stop("No saved data found. Run Script 02 first.")
#       01_Feb2.R        -> readRDS error on first_stage_pm25.rds
#       01_Feb2_3.R      -> falls through to a partial fallback branch
#
# THE FIX
#   Read the "_2" objects once and re-save them under every name the
#   downstream scripts expect. No existing script is edited.
#
# RUN THIS IMMEDIATELY AFTER 0202_first.R AND BEFORE ANYTHING ELSE.
################################################################################

cat("################################################################\n")
cat("# 00_Bridge_RDS.R - aligning RDS filenames\n")
cat("################################################################\n\n")

RDS <- "Output/RDS_files"
if (!dir.exists(RDS)) stop("Output/RDS_files does not exist. Run 0202_first.R first.")

#-------------------------------------------------------------------------------
# 1. Load the objects 0202_first.R actually wrote
#-------------------------------------------------------------------------------

need <- c(
  merged_daily_pm2.5 = "merged_daily_pm25_2.rds",
  merged_daily_pm10  = "merged_daily_pm10_2.rds",
  results_pm25       = "results_pm25_2.rds",
  results_pm10       = "results_pm10_2.rds",
  fs_pm25            = "first_stage_pm25_2.rds",
  fs_pm10            = "first_stage_pm10_2.rds",
  ols_pm25           = "ols_pm25_2.rds",
  ols_pm10           = "ols_pm10_2.rds",
  dayfe_pm25         = "dayfe_pm25_2.rds",
  dayfe_pm10         = "dayfe_pm10_2.rds"
)

missing <- character(0)
for (nm in names(need)) {
  f <- file.path(RDS, need[[nm]])
  if (file.exists(f)) {
    assign(nm, readRDS(f))
    cat("  loaded ", need[[nm]], "\n", sep = "")
  } else {
    missing <- c(missing, need[[nm]])
  }
}

if (length(missing)) {
  cat("\n  WARNING - not found (0202_first.R may not have completed):\n")
  for (m in missing) cat("    ", m, "\n")
  cat("\n")
}

#-------------------------------------------------------------------------------
# 2. Re-save under the unsuffixed names the 01_Feb2 scripts read
#-------------------------------------------------------------------------------

alias <- list(
  merged_daily_pm2.5 = c("merged_daily_pm25.rds", "data_pm25.rds"),
  merged_daily_pm10  = c("merged_daily_pm10.rds", "data_pm10.rds"),
  results_pm25       = "results_pm25.rds",
  results_pm10       = "results_pm10.rds",
  fs_pm25            = c("first_stage_pm25.rds", "fs_pm25.rds"),
  fs_pm10            = c("first_stage_pm10.rds", "fs_pm10.rds"),
  ols_pm25           = "ols_pm25.rds",
  ols_pm10           = "ols_pm10.rds",
  dayfe_pm25         = "dayfe_pm25.rds",
  dayfe_pm10         = "dayfe_pm10.rds"
)

cat("\n  Writing aliases:\n")
for (nm in names(alias)) {
  if (!exists(nm)) next
  for (target in alias[[nm]]) {
    saveRDS(get(nm), file.path(RDS, target))
    cat("    ", target, "\n", sep = "")
  }
}

#-------------------------------------------------------------------------------
# 3. Write the .RData bundle that 01_Feb2_2.R and 01_Feb2_4_PM10.R require
#-------------------------------------------------------------------------------

bundle <- intersect(
  c("merged_daily_pm2.5", "merged_daily_pm10", "results_pm25", "results_pm10",
    "fs_pm25", "fs_pm10", "ols_pm25", "ols_pm10", "dayfe_pm25", "dayfe_pm10"),
  ls()
)

save(list = bundle, file = file.path(RDS, "script02_all_results.RData"))
cat("\n  Wrote script02_all_results.RData containing:\n    ",
    paste(bundle, collapse = ", "), "\n", sep = "")

#-------------------------------------------------------------------------------
# 4. Sanity check on the SE field name
#-------------------------------------------------------------------------------
# 01_Feb2_3.R probes results_pm25$se_twoway, then $twoway_se, then
# $bootstrap_se, then $cf_se, and silently proceeds with whichever it finds
# first. If the field it lands on is not the two-way clustered SE, every
# downstream confidence interval is built on the wrong variance.

if (exists("results_pm25")) {
  cat("\n  Fields available in results_pm25:\n    ",
      paste(names(results_pm25), collapse = ", "), "\n", sep = "")
  cat("\n  CHECK: confirm the two-way clustered SE is stored under 'se_twoway'.\n")
  cat("  If it is stored under a different name, 01_Feb2_3.R will pick up a\n")
  cat("  bootstrap or unclustered SE instead without warning.\n")
}

cat("\nBridge complete. Downstream scripts can now find their inputs.\n")
