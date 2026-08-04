################################################################################
# 92_NegativeBinomial.R
#
# GAP 2 of 5  ->  MANUSCRIPT TABLE 11, ROW "Negative binomial (alternative GLM)"
#
# Target values to reproduce:  PM2.5  beta = 0.0152, SE = 0.0039
#
# Why this script exists: glm.nb / MASS / "negative binomial" appears in ZERO
# of the 24 scripts. The Table 11 row has no code behind it.
#
# What it does:
#   1. Tests for overdispersion in the Poisson second stage (the thing that
#      motivates a negative binomial in the first place)
#   2. Re-estimates the control-function second stage as NB2 via MASS::glm.nb
#   3. Reports the Poisson benchmark alongside, so the row is interpretable
#
# RUNTIME WARNING: glm.nb with ~60 fixed-effect dummies on 55,719 rows takes
# roughly 3-10 minutes per pollutant. This is expected. Do not kill it.
#
# Output: Output/Gap2_negative_binomial.csv
################################################################################

source("90_Common_Setup.R")

cat("################################################################\n")
cat("# GAP 2: NEGATIVE BINOMIAL SECOND STAGE\n")
cat("################################################################\n\n")

merged_daily <- prep_analysis_data()

results <- list()

for (poll in c("PM2.5AVG", "PM10AVG")) {
  
  cat("---------------------------------------------------------------\n")
  cat(poll, "\n")
  cat("---------------------------------------------------------------\n")
  
  d <- complete_sample(merged_daily, poll)
  
  #---------------------------------------------------------------------------
  # 1. POISSON BENCHMARK
  #---------------------------------------------------------------------------
  r_pois <- run_cf_iv(d, poll, family = "poisson", cluster = "twoway",
                      label = paste(poll, "| Poisson (benchmark)"))
  
  #---------------------------------------------------------------------------
  # 2. OVERDISPERSION TEST
  #    Justifies the NB. If dispersion ~= 1, say so explicitly in the note -
  #    it strengthens the paper rather than weakening it.
  #---------------------------------------------------------------------------
  if (isTRUE(r_pois$ok)) {
    ss   <- r_pois$second_stage
    pr   <- residuals(ss, type = "pearson")
    disp <- sum(pr^2) / df.residual(ss)
    
    cat(sprintf("  Overdispersion (Pearson chi2 / df) : %.3f\n", disp))
    cat(sprintf("  Outcome mean = %.3f, variance = %.3f\n",
                mean(d$deaths), var(d$deaths)))
    if (disp > 1.2) {
      cat("  => Overdispersed. Negative binomial is the appropriate check.\n\n")
    } else {
      cat("  => Close to equidispersion. Poisson is already adequate;\n")
      cat("     the NB row then serves as reassurance, not correction.\n\n")
    }
  } else {
    disp <- NA_real_
  }
  
  #---------------------------------------------------------------------------
  # 3. NEGATIVE BINOMIAL SECOND STAGE
  #---------------------------------------------------------------------------
  cat("  Estimating NB2 (this is slow - several minutes)...\n")
  t0 <- Sys.time()
  
  r_nb <- run_cf_iv(d, poll, family = "negbin", cluster = "twoway",
                    label = paste(poll, "| Negative binomial"))
  
  cat(sprintf("  Elapsed: %.1f minutes\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  if (isTRUE(r_nb$ok)) {
    cat(sprintf("  theta = %.3f (SE %.3f)   [theta -> Inf means NB collapses to Poisson]\n\n",
                r_nb$theta, r_nb$theta_se))
  }
  
  results[[paste0(poll, "_poisson")]] <- r_pois
  results[[paste0(poll, "_negbin")]]  <- r_nb
  
  if (isTRUE(r_pois$ok) && isTRUE(r_nb$ok)) {
    cat(sprintf("  => Poisson vs NB coefficient difference: %+.2f%%\n\n",
                100 * (r_nb$coefficient - r_pois$coefficient) / r_pois$coefficient))
  }
}

#===============================================================================
# ASSEMBLE TABLE 11 ROW
#===============================================================================

tab <- do.call(rbind, lapply(results, cf_row))
tab$Theta <- sapply(results, function(x) if (isTRUE(x$ok)) x$theta else NA)
rownames(tab) <- NULL

cat("\n===============================================================\n")
cat("TABLE 11 INPUT - NEGATIVE BINOMIAL\n")
cat("===============================================================\n")
print(tab %>%
        mutate(across(c(Coefficient, SE), ~ round(.x, 5)),
               across(c(Effect_per10, Theta), ~ round(.x, 2)),
               P_value = round(P_value, 4)) %>%
        dplyr::select(Label, Coefficient, SE, P_value, Effect_per10, Theta, N))

write.csv(tab, "Output/Gap2_negative_binomial.csv", row.names = FALSE)
cat("\nSaved: Output/Gap2_negative_binomial.csv\n")

r <- results[["PM2.5AVG_negbin"]]
if (isTRUE(r$ok)) {
  cat("\n--- Reconciliation with Table 11 ---\n")
  cat(sprintf("  Manuscript : beta = 0.0152, SE = 0.0039\n"))
  cat(sprintf("  This run   : beta = %.4f, SE = %.4f\n", r$coefficient, r$se))
  cat("\n  NOTE: the F-stat column for this row in Table 11 currently shows\n")
  cat("  5,658, i.e. the main-specification first stage. That is correct -\n")
  cat("  the NB changes only the second stage - but the table note should\n")
  cat("  say so, or a referee will read it as a separately estimated F.\n")
  cat("\n  CAVEAT for the note: the two-way clustered SE treats theta as\n")
  cat("  fixed at its estimated value. This is standard, but state it.\n")
}

cat("\nGap 2 complete.\n")


