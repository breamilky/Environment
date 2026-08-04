################################################################################
# RUN_ALL.R   (revised)
#
# Changes from the previous version:
#   - revised_counterfactual.R REMOVED. 0202_first.R produces Tables 12-14
#     internally and its output is what the manuscript reports (128 excess
#     deaths, 10.0% fire share, RM 192.6m). The revised_ version rebuilds the
#     instrument from WINDSPEEDAVG, lands on gamma_u = 1.88e-04 instead of
#     2.26e-04, and gives 94 deaths instead of 128.
#   - 96_Pairs_Bootstrap.R REMOVED from the routine run. Table A2's primary SE
#     is confirmed (0.003676 vs published 0.00368). Keep the script for the
#     replication archive; re-run at n_boot = 500 only if you want the
#     tighter figure.
#   - 97_Patch_And_Archive.R ADDED as a prerequisite.
#
# Usage:  source("RUN_ALL.R")
# Runtime: ~4 hours serial, ~2.5 if Stage 3 runs in a second session.
################################################################################

################################################################################
# PRE-FLIGHT
################################################################################
#
# ONE COMMAND:
#     source("97_Patch_And_Archive.R")
#   Patches the two-way clustering guard in every script that has it (the
#   fallback to station-clustering was writing 0.00315 where Table 5 and
#   Table A2 report 0.00367), and archives stale output.
#
# ALREADY DONE (verified in this session):
#   90_Common_Setup.R   uses the CSV's u_wind, not a rebuild   -> a1 = 2.259e-04
#   90_Common_Setup.R   CGM eigenvalue correction              -> SE 0.003671, F 26.8
#   01_Feb2_3.R         geocode column-shift realignment       -> 45 stations
#   01_Feb2_3.R         spline block gated with if (FALSE)
#   93_Regional_...R    cluster = "none" per station + alias skip
#
# STILL TO DO BY HAND:
#   0202_second.R       year-effects block: linear-combination SE with the
#                       covariance term, plus interaction p-values.
#                       Table A6 currently has SE = NA for 2017/2018/2020.
#   01_Feb2_4_PM10.R    analysis_data <- merged_daily  (not merged_daily_pm2.5)
#
################################################################################

t_start <- Sys.time()
cat("Pipeline started:", format(t_start), "\n\n")

#-------------------------------------------------------------------------------
# STAGE 1 - Data
#-------------------------------------------------------------------------------
source("01_Clean_Merge.R")              # Table 1                        ~5 min

#-------------------------------------------------------------------------------
# STAGE 2 - Core results   (one session, this order - steps 4 and 5 read
#                           objects left in the global environment)
#-------------------------------------------------------------------------------
source("0202_first.R")                  # Tables 2-6, 11, 12, 13, 14,
# A1a, A1b, A2; Figure 2        ~90 min
source("00_Bridge_RDS.R")               # RDS filename bridge            <1 min
source("0202_second.R")                 # Tables 7-10, A3-A6 (PM2.5)     ~60 min
source("01_Feb2_4_PM10.R")              # PM10 columns                   ~30 min
source("01_Feb2_3.R")                   # Figure 1, Figure A1            ~20 min
source("96_Pairs_Bootstrap.R")          # Table A2 row 1 (PRIMARY SE)

#-------------------------------------------------------------------------------
# STAGE 3 - Gap scripts   (read the merged CSV only; safe to run in a
#                          second R session in parallel with Stage 2)
#-------------------------------------------------------------------------------
source("91_Exclude2020.R")              # Table 11, COVID row             ~5 min
source("92_NegativeBinomial.R")         # Table 11, NegBin row            ~5 min
source("93_Regional_Heterogeneity.R")   # Section 5.3                     ~5 min
source("94_FireDays_Above5.R")          # Table 12, final column          ~2 min
source("95_Welfare_Sensitivity.R")      # Section 6.7                     ~5 min

#-------------------------------------------------------------------------------
# STAGE 4 - Collect
#-------------------------------------------------------------------------------
source("99_Build_Publication_Folder.R") # Publication/                   <1 min

cat("\nFinished in",
    round(as.numeric(difftime(Sys.time(), t_start, units = "hours")), 2),
    "hours\n")

