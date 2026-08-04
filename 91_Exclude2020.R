
################################################################################
# 91_Exclude2020.R
#
# GAP 1 of 5  ->  MANUSCRIPT TABLE 11, ROW "Excluding 2020 (COVID period)"
#
# Target values to reproduce:  PM2.5  beta = 0.0161, SE = 0.0040, F = 698
#
# Why this script exists: no script in the existing codebase filters 2020
# anywhere. The claim in Section 7 ("Our results are robust to excluding
# 2020") is currently unsupported by any code.
#
# Also produces the PM10 equivalent, which the manuscript does not report but
# a referee will almost certainly ask for.
#
# Output: Output/Gap1_exclude2020.csv
################################################################################

source("90_Common_Setup.R")

cat("################################################################\n")
cat("# GAP 1: ROBUSTNESS TO EXCLUDING 2020 (COVID PERIOD)\n")
cat("################################################################\n\n")

merged_daily <- prep_analysis_data()

results <- list()

for (poll in c("PM2.5AVG", "PM10AVG")) {
  
  cat("---------------------------------------------------------------\n")
  cat(poll, "\n")
  cat("---------------------------------------------------------------\n")
  
  d_full <- complete_sample(merged_daily, poll)
  d_excl <- d_full %>% filter(year != 2020)
  
  cat("  Full sample      :", format(nrow(d_full), big.mark = ","), "obs,",
      format(sum(d_full$deaths), big.mark = ","), "deaths\n")
  cat("  Excluding 2020   :", format(nrow(d_excl), big.mark = ","), "obs,",
      format(sum(d_excl$deaths), big.mark = ","), "deaths",
      sprintf("(%.1f%% of sample dropped)\n\n",
              100 * (1 - nrow(d_excl) / nrow(d_full))))
  
  # Benchmark: full sample, so the two rows sit side by side in the table
  r_full <- run_cf_iv(d_full, poll, cluster = "twoway",
                      label = paste(poll, "| full sample (benchmark)"))
  
  # The row the manuscript actually reports
  r_excl <- run_cf_iv(d_excl, poll, cluster = "twoway",
                      label = paste(poll, "| excluding 2020"))
  
  results[[paste0(poll, "_full")]] <- r_full
  results[[paste0(poll, "_excl")]] <- r_excl
  
  if (isTRUE(r_full$ok) && isTRUE(r_excl$ok)) {
    delta <- 100 * (r_excl$coefficient - r_full$coefficient) / r_full$coefficient
    cat(sprintf("  => Coefficient change from dropping 2020: %+.1f%%\n\n", delta))
  }
}

#===============================================================================
# ASSEMBLE TABLE 11 ROWS
#===============================================================================

tab <- do.call(rbind, lapply(results, cf_row))
rownames(tab) <- NULL

cat("\n===============================================================\n")
cat("TABLE 11 INPUT - COVID ROBUSTNESS\n")
cat("===============================================================\n")
print(tab %>%
        mutate(across(c(Coefficient, SE), ~ round(.x, 5)),
               across(c(Effect_per10, F_unclustered, F_clustered), ~ round(.x, 1)),
               P_value = round(P_value, 4)) %>%
        dplyr::select(Label, Coefficient, SE, P_value, Effect_per10,
                      F_unclustered, F_clustered, N))

write.csv(tab, "Output/Gap1_exclude2020.csv", row.names = FALSE)
cat("\nSaved: Output/Gap1_exclude2020.csv\n")

#===============================================================================
# CHECK AGAINST THE MANUSCRIPT
#===============================================================================

r <- results[["PM2.5AVG_excl"]]
if (isTRUE(r$ok)) {
  cat("\n--- Reconciliation with Table 11 ---\n")
  cat(sprintf("  Manuscript : beta = 0.0161, SE = 0.0040, F = 698\n"))
  cat(sprintf("  This run   : beta = %.4f, SE = %.4f, F = %.0f (unclustered)\n",
              r$coefficient, r$se, r$f_unclustered))
  cat(sprintf("                                     F = %.1f (two-way clustered)\n",
              r$f_clustered))
}

cat("\nGap 1 complete.\n")




