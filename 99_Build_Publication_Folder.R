################################################################################
# 99_Build_Publication_Folder.R
#
# Collects every table and figure the manuscript reports into Publication/,
# renamed to match the manuscript's own numbering, and reports what is
# missing.
#
# The pipeline scatters output across Output/, Output/Robustness/Tables/,
# Output/PM10_Robustness/Tables/, Output/DoseResponse/ and Output/Figures/,
# with filenames that do not indicate which table they feed. This script does
# not re-estimate anything - it copies finished output into one place.
#
# RUN LAST, after RUN_ALL.R completes.
#
# Usage:  source("99_Build_Publication_Folder.R")
################################################################################

cat("################################################################\n")
cat("# BUILDING Publication/\n")
cat("################################################################\n\n")

PUB <- "Publication"

for (sub in c("Tables", "Figures", "Supporting")) {
  dir.create(file.path(PUB, sub), recursive = TRUE, showWarnings = FALSE)
}

#===============================================================================
# MANIFEST
#===============================================================================
# dest      = filename inside Publication/
# src       = candidate source paths, first match wins
# manuscript= what it feeds
# critical  = TRUE if the manuscript cannot be submitted without it

manifest <- list(
  
  #--- MAIN TABLES ------------------------------------------------------------
  list(dest = "Tables/Table_01_Data_Structure.csv",
       src  = c("Output/tables/Table1_Panel_Structure_2.csv"),
       ms   = "Table 1: Data structure (45 pairs, 55,719 obs, 1,277 days)",
       critical = TRUE),
  
  list(dest = "Tables/Table_02_Summary_Statistics.csv",
       src  = c("Output/descriptive_statistics_2.csv",
                "Output/Tables/descriptive_statistics_2.csv"),
       ms   = "Table 2: Summary statistics",
       critical = TRUE),
  
  list(dest = "Tables/Table_03_First_Stage.csv",
       src  = c("Output/first_stage_table_2.csv",
                "Output/Tables/first_stage_table_2.csv"),
       ms   = "Table 3: First-stage regression",
       critical = TRUE),
  
  list(dest = "Tables/Table_04_Reduced_Form.csv",
       src  = c("Output/tables/Table6A_Reduced_Form_Coefficients_2.csv"),
       ms   = "Table 4: Reduced-form (instruments -> mortality)",
       critical = TRUE),
  
  list(dest = "Tables/Table_04b_Reduced_Form_Joint_Test.csv",
       src  = c("Output/tables/Table6B_Reduced_Form_Joint_Test_2.csv"),
       ms   = "Table 4: joint significance of instruments",
       critical = TRUE),
  
  list(dest = "Tables/Table_05_Main_Results.csv",
       src  = c("Output/main_results_table_2.csv",
                "Output/Tables/main_results_table_2.csv"),
       ms   = "Table 5: OLS vs IV main results",
       critical = TRUE),
  
  list(dest = "Tables/Table_06_Variance_Decomposition.csv",
       src  = c("Output/tables/Table7A_Variance_Decomposition_2.csv"),
       ms   = "Table 6: Variance decomposition",
       critical = TRUE),
  
  list(dest = "Tables/Table_07_Placebo.csv",
       src  = c("Output/Robustness/Tables/Placebo_Outcomes.csv",
                "Output/Placebo_Outcomes.csv"),
       ms   = "Table 7: Placebo outcome tests",
       critical = TRUE),
  
  list(dest = "Tables/Table_08_Season_Heterogeneity.csv",
       src  = c("Output/Robustness/Tables/Heterogeneity_Season.csv",
                "Output/Heterogeneity_Season.csv"),
       ms   = "Table 8: Seasonal heterogeneity",
       critical = TRUE),
  
  list(dest = "Tables/Table_09_Dose_Response.csv",
       src  = c("Output/DoseResponse/dose_response_bins.csv",
                "Output/dose_response_bins.csv"),
       ms   = "Table 9: Dose-response by PM2.5 bin",
       critical = TRUE),
  
  list(dest = "Tables/Table_10_Cumulative_Effects.csv",
       src  = c("Output/Robustness/Tables/Cumulative_Effects.csv",
                "Output/Cumulative_Effects.csv"),
       ms   = "Table 10: Cumulative exposure effects",
       critical = TRUE),
  
  list(dest = "Tables/Table_11a_FE_Robustness.csv",
       src  = c("Output/fe_comparison_table_2.csv",
                "Output/Tables/fe_comparison_table_2.csv"),
       ms   = "Table 11: rows 1-2, 5-6 (FE specifications, day FE)",
       critical = TRUE),
  
  list(dest = "Tables/Table_11b_Exclude_2020.csv",
       src  = c("Output/Gap1_exclude2020.csv"),
       ms   = "Table 11: row 3 (excluding 2020, COVID)",
       critical = TRUE),
  
  list(dest = "Tables/Table_11c_Negative_Binomial.csv",
       src  = c("Output/Gap2_negative_binomial.csv"),
       ms   = "Table 11: row 4 (negative binomial)",
       critical = TRUE),
  
  list(dest = "Tables/Table_12_Fire_Attributable.csv",
       src  = c("Output/Gap4_table12_fire_attributable.csv"),
       ms   = "Table 12: fire-attributable pollution incl. days > 5 ug/m3",
       critical = TRUE),
  
  list(dest = "Tables/Table_13_Counterfactual_By_Year.csv",
       src  = c("Output/Gap5_table13_rebuilt.csv",
                "Output/counterfactual_pm25_by_year.csv",
                "Output/counterfactual_by_year.csv"),
       ms   = "Table 13: excess deaths by year (per-year denominator)",
       critical = TRUE),
  
  list(dest = "Tables/Table_14_Welfare.csv",
       src  = c("Output/Gap5_welfare_sensitivity.csv",
                "Output/counterfactual_summary_corrected.csv"),
       ms   = "Table 14: welfare estimates + VSL sensitivity",
       critical = TRUE),
  
  #--- APPENDIX TABLES --------------------------------------------------------
  list(dest = "Tables/Table_A1a_FE_Sensitivity.csv",
       src  = c("Output/fe_comparison_table_2.csv"),
       ms   = "Table A1a: sensitivity to FE specification",
       critical = TRUE),
  
  list(dest = "Tables/Table_A1b_Power_Analysis.csv",
       src  = c("Output/power_analysis_2.csv",
                "Output/power_analysis.csv",
                "Output/Tables/power_analysis_2.csv"),
       ms   = "Table A1b: power analysis / MDE",
       critical = TRUE),
  
  list(dest = "Tables/Table_A2_Bootstrap_Sensitivity.csv",
       src  = c("Output/TableA2_pairs_bootstrap.csv"),
       ms   = "Table A2: bootstrap sensitivity (PRIMARY inference)",
       critical = TRUE),
  
  list(dest = "Tables/Table_A3_Distributed_Lag.csv",
       src  = c("Output/Robustness/Tables/Distributed_Lag_Model.csv",
                "Output/Distributed_Lag_Model.csv"),
       ms   = "Table A3: distributed lag structure",
       critical = TRUE),
  
  list(dest = "Tables/Table_A4_Alternative_IV.csv",
       src  = c("Output/Robustness/Tables/Alternative_IV_Specifications.csv",
                "Output/Alternative_IV_Specifications.csv"),
       ms   = "Table A4: alternative instrument constructions",
       critical = TRUE),
  
  list(dest = "Tables/Table_A5_Outlier_Sensitivity.csv",
       src  = c("Output/Robustness/Tables/Sensitivity_Outliers.csv",
                "Output/Sensitivity_Outliers.csv"),
       ms   = "Table A5: outlier / winsorising sensitivity",
       critical = TRUE),
  
  list(dest = "Tables/Table_A6_Year_Heterogeneity_PM25.csv",
       src  = c("Output/Robustness/Tables/Heterogeneity_Year.csv",
                "Output/Heterogeneity_Year.csv"),
       ms   = "Table A6: year heterogeneity, PM2.5 column",
       critical = TRUE),
  
  list(dest = "Tables/Table_A6_Year_Heterogeneity_PM10.csv",
       src  = c("Output/PM10_Robustness/Tables/PM10_Heterogeneity_Year.csv"),
       ms   = "Table A6: year heterogeneity, PM10 column",
       critical = TRUE),
  
  #--- SECTION 5.3 (no table number, reported in text) ------------------------
  list(dest = "Tables/Section_5-3_Regional_Heterogeneity.csv",
       src  = c("Output/Gap3_regional_heterogeneity.csv"),
       ms   = "Section 5.3: regional heterogeneity, 45 stations",
       critical = TRUE),
  
  #--- FIGURES ----------------------------------------------------------------
  list(dest = "Figures/Figure_01_Within_Day_Variation.png",
       src  = c("Output/Maps/identifying_variation_map.png"),
       ms   = "Figure 1: within-day identifying variation, 11 Sep 2019",
       critical = TRUE),
  
  list(dest = "Figures/Figure_A1_Event_Study.png",
       src  = c("Output/event_study_improved.png"),
       ms   = "Figure A1: event study around severe pollution episodes",
       critical = TRUE),
  
  #--- SUPPORTING (not in the manuscript, but referee-facing) ------------------
  list(dest = "Supporting/Event_Study_Coefficients.csv",
       src  = c("Output/event_study_coefficients_improved.csv",
                "Output/DoseResponse/event_study_coefficients.csv"),
       ms   = "Figure A1 underlying coefficients",
       critical = FALSE),
  
  list(dest = "Supporting/Event_Study_Threshold_Sensitivity.csv",
       src  = c("Output/event_study_threshold_sensitivity.csv"),
       ms   = "Section 4.5: 85th/90th/95th percentile sensitivity",
       critical = FALSE),
  
  list(dest = "Supporting/Regional_Forest_PM10.png",
       src  = c("Output/Gap3_regional_forest_PM10.png"),
       ms   = "Section 5.3 forest plot, PM10",
       critical = FALSE),
  
  list(dest = "Supporting/Regional_Forest_PM25.png",
       src  = c("Output/Gap3_regional_forest_PM25.png"),
       ms   = "Section 5.3 forest plot, PM2.5",
       critical = FALSE),
  
  list(dest = "Supporting/Bootstrap_Draws.rds",
       src  = c("Output/TableA2_bootstrap_draws.rds"),
       ms   = "Table A2 replicate distribution, for diagnostics",
       critical = FALSE),
  
  list(dest = "Supporting/Dose_Response_Spline.csv",
       src  = c("Output/DoseResponse/dose_response_spline.csv"),
       ms   = "Nonlinearity check, not reported",
       critical = FALSE),
  
  list(dest = "Supporting/Lead_Placebo.csv",
       src  = c("Output/Robustness/Tables/lead_placebo_summary_improved.csv",
                "Output/lead_placebo_summary_improved.csv"),
       ms   = "Lead placebo / pre-trend check",
       critical = FALSE),
  
  list(dest = "Supporting/Threshold_Tests.csv",
       src  = c("Output/threshold_test_results.csv"),
       ms   = "Threshold / nonlinearity tests",
       critical = FALSE)
)

#===============================================================================
# COPY
#===============================================================================

log <- data.frame(Destination = character(), Manuscript = character(),
                  Status = character(), Source = character(),
                  Critical = logical(), stringsAsFactors = FALSE)

for (item in manifest) {
  
  hit <- item$src[file.exists(item$src)][1]
  
  if (!is.na(hit)) {
    file.copy(hit, file.path(PUB, item$dest), overwrite = TRUE)
    status <- "OK"
    cat(sprintf("  OK       %-46s <- %s\n", item$dest, hit))
  } else {
    status <- if (isTRUE(item$critical)) "MISSING (CRITICAL)" else "missing"
    hit    <- NA_character_
    cat(sprintf("  %-8s %-46s\n",
                if (isTRUE(item$critical)) "MISSING!" else "skip", item$dest))
  }
  
  log <- rbind(log, data.frame(
    Destination = item$dest,
    Manuscript  = item$ms,
    Status      = status,
    Source      = ifelse(is.na(hit), "", hit),
    Critical    = isTRUE(item$critical),
    stringsAsFactors = FALSE
  ))
}

#===============================================================================
# MANIFEST FILE
#===============================================================================

write.csv(log, file.path(PUB, "MANIFEST.csv"), row.names = FALSE)

n_ok       <- sum(log$Status == "OK")
n_crit_bad <- sum(log$Status == "MISSING (CRITICAL)")

cat("\n===============================================================\n")
cat(sprintf("Collected %d of %d items into %s/\n", n_ok, nrow(log), PUB))
cat(sprintf("Critical items missing: %d\n", n_crit_bad))
cat("===============================================================\n")

if (n_crit_bad > 0) {
  cat("\nMISSING CRITICAL OUTPUT - the manuscript cannot be assembled:\n\n")
  bad <- log[log$Status == "MISSING (CRITICAL)", ]
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  %s\n      %s\n", bad$Destination[i], bad$Manuscript[i]))
  }
  cat("\nIf a file exists under a different name, add its path to the\n")
  cat("`src` vector for that entry and re-run. The candidate paths are\n")
  cat("guesses for the scripts whose exact output names I could not verify\n")
  cat("(Tables 1, 4, 6, A1b and both figures).\n")
}

cat("\nSaved: ", PUB, "/MANIFEST.csv\n", sep = "")
cat("\nNOTE: Table 11 and Table A1a both draw on fe_comparison_table_2.csv;\n")
cat("it is copied twice under both names so each table has its own file.\n")