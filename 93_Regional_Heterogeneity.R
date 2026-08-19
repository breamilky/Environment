################################################################################
# 93_Regional_Heterogeneity.R
#
# GAP 3 of 5  ->  MANUSCRIPT SECTION 5.3 "Regional Heterogeneity"
#
# Target values to reproduce:
#   PM10 coefficients 0.0506 (Alor Setar) and 0.0508 (Sungai Petani)
#   vs sample-wide 0.0141; Pegoh Ipoh and Batu Muda "42-161% above average"
#
# Why this script exists: the ONLY regional code you have is 03_Region_Lag.R,
# which uses the superseded specification - linear ivreg, three instruments
# (total_frp + non-station-specific FRP_u_wind + FRP_v_wind), no station FE,
# Newey-West SEs. Section 5.3 claims IV-CF estimates. Those do not exist yet.
#
# This is the largest of the five gaps: 45 separate two-stage estimations.
#
# NOTE ON WITHIN-STATION SPECIFICATION: with a single station there is no
# cross-sectional dimension, so station FE drop out and clustering is by date
# only. Everything else matches the main spec.
#
# Runtime: roughly 2-5 minutes for all 45 stations x 2 pollutants.
#
# Outputs: Output/Gap3_regional_heterogeneity.csv
#          Output/Gap3_regional_forest_PM10.png
#          Output/Gap3_regional_forest_PM25.png
################################################################################

source("90_Common_Setup.R")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

cat("################################################################\n")
cat("# GAP 3: REGIONAL HETEROGENEITY (45 STATIONS, IV-CF)\n")
cat("################################################################\n\n")

merged_daily <- prep_analysis_data()

# Minimum observations before a station is attempted at all
MIN_OBS    <- 200
MIN_DEATHS <- 30

regional_all <- list()

for (poll in c("PM10AVG", "PM2.5AVG")) {
  
  cat("===============================================================\n")
  cat("REGIONAL ESTIMATES:", poll, "\n")
  cat("===============================================================\n\n")
  
  d_all <- complete_sample(merged_daily, poll)
  
  #---------------------------------------------------------------------------
  # National benchmark (the number Section 5.3 compares each station against)
  #---------------------------------------------------------------------------
  nat <- run_cf_iv(d_all, poll, cluster = "twoway",
                   label = paste(poll, "| NATIONAL (sample-wide)"))
  beta_nat <- nat$coefficient
  
  #---------------------------------------------------------------------------
  # Station-by-station
  #---------------------------------------------------------------------------
  stations <- sort(unique(as.character(d_all$station_id)))
  cat("Estimating", length(stations), "stations...\n\n")
  
  rows <- list()
  
  for (st in stations) {
    
    d_st <- d_all[as.character(d_all$station_id) == st, , drop = FALSE]
    
    if (nrow(d_st) < MIN_OBS || sum(d_st$deaths) < MIN_DEATHS) {
      cat(sprintf("  %-32s SKIPPED (n=%d, deaths=%d)\n",
                  st, nrow(d_st), sum(d_st$deaths)))
      next
    }
    
    # Single station: no station FE, cluster by date only
    r <- run_cf_iv(d_st, poll,
                   use_station_fe = FALSE,
                   cluster        = "none",   
                   label          = st,          # see note below
                   quiet          = TRUE)
    
    # For a single station, station_id is constant, so "station" clustering
    # degenerates to a single cluster. Re-do inference clustering on DATE,
    # which is the correct unit here.
    if (isTRUE(r$ok)) {
      d_fit <- complete_sample(d_st, poll)
      d_fit$cf_resid <- residuals(r$first_stage)
      vc <- sandwich::vcovCL(r$second_stage,
                             cluster = as.character(d_fit$date), type = "HC1")
      b  <- coef(r$second_stage)[poll]
      se <- sqrt(vc[poll, poll])
      
      r$coefficient  <- unname(b)
      r$se           <- unname(se)
      r$p_value      <- unname(2 * pnorm(-abs(b / se)))
      r$ci_lo        <- unname(b - 1.96 * se)
      r$ci_hi        <- unname(b + 1.96 * se)
      r$effect_pct10 <- (exp(unname(b) * 10) - 1) * 100
    }
    
    if (!isTRUE(r$ok)) {
      cat(sprintf("  %-32s FAILED: %s\n", st, r$reason))
      next
    }
    
    # Two ways of expressing the comparison to national. Section 5.3 currently
    # conflates them - see the diagnostic printed at the end of this script.
    ratio_of    <- r$coefficient / beta_nat            # "X% OF national"
    pct_above   <- 100 * (r$coefficient - beta_nat) / beta_nat  # "X% ABOVE"
    
    rows[[st]] <- data.frame(
      Pollutant      = poll,
      Station        = st,
      N              = r$n,
      Deaths         = r$deaths,
      Coefficient    = r$coefficient,
      SE             = r$se,
      P_value        = r$p_value,
      CI_Lower       = r$ci_lo,
      CI_Upper       = r$ci_hi,
      Effect_per10   = r$effect_pct10,
      First_stage_F  = r$f_unclustered,
      National_beta  = beta_nat,
      Pct_OF_national    = 100 * ratio_of,
      Pct_ABOVE_national = pct_above,
      Significant    = r$ci_lo > 0 | r$ci_hi < 0,
      Weak_IV        = r$f_unclustered < 10,
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("  %-32s beta = %8.5f (SE %7.5f)  p=%.3f  F=%6.0f  %s\n",
                st, r$coefficient, r$se, r$p_value, r$f_unclustered,
                ifelse(r$f_unclustered < 10, "<- WEAK IV", "")))
  }
  
  reg <- do.call(rbind, rows)

  # An empty result should report itself, not crash at the arrange() below.
  if (is.null(reg) || nrow(reg) == 0) {
    cat("\n  No stations estimated for", poll, "\n")
    cat("  Thresholds: MIN_OBS =", MIN_OBS, ", MIN_DEATHS =", MIN_DEATHS, "\n")
    cat("  On a subsample, lower these at lines 42-43 (try 30 and 3).\n\n")
    next
  }
  
  rownames(reg) <- NULL
  regional_all[[poll]] <- reg
  
  #---------------------------------------------------------------------------
  # Summary for this pollutant
  #---------------------------------------------------------------------------
  cat("\n  ---- Summary ----\n")
  cat(sprintf("  Stations estimated      : %d of %d\n",
              nrow(reg), length(stations)))
  cat(sprintf("  Statistically significant: %d\n", sum(reg$Significant)))
  cat(sprintf("  Weak first stage (F<10) : %d\n", sum(reg$Weak_IV)))
  cat(sprintf("  National benchmark      : %.5f\n", beta_nat))
  cat(sprintf("  Range across stations   : %.5f to %.5f\n",
              min(reg$Coefficient), max(reg$Coefficient)))
  
  cat("\n  Top 5 stations by effect size:\n")
  print(reg %>%
          arrange(desc(Coefficient)) %>%
          head(5) %>%
          mutate(Coefficient = round(Coefficient, 5),
                 SE = round(SE, 5),
                 Pct_OF_national = round(Pct_OF_national, 0),
                 Pct_ABOVE_national = round(Pct_ABOVE_national, 0)) %>%
          dplyr::select(Station, Coefficient, SE, P_value,
                        Pct_OF_national, Pct_ABOVE_national, First_stage_F))
  
  #---------------------------------------------------------------------------
  # Forest plot
  #---------------------------------------------------------------------------
  p <- ggplot(reg, aes(x = reorder(Station, Coefficient), y = Coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = beta_nat, linetype = "solid", colour = "firebrick") +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                  width = 0.25, colour = "grey30") +
    geom_point(aes(colour = Significant), size = 2) +
    scale_colour_manual(values = c(`TRUE` = "steelblue4", `FALSE` = "grey60"),
                        name = "95% CI excludes 0") +
    coord_flip() +
    labs(
      title    = paste("Regional heterogeneity in", poll, "mortality effects"),
      subtitle = sprintf("IV-CF estimates by monitoring station. Red line = national estimate (%.4f)",
                         beta_nat),
      x = NULL, y = "Coefficient (log mortality per ug/m3)"
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "bottom")
  
  fn <- paste0("Output/Gap3_regional_forest_",
               ifelse(poll == "PM10AVG", "PM10", "PM25"), ".png")
  ggsave(fn, p, width = 8, height = 10, dpi = 300)
  cat("\n  Saved plot:", fn, "\n\n")
}

#===============================================================================
# COMBINED OUTPUT
#===============================================================================

out <- do.call(rbind, regional_all)

if (is.null(out) || nrow(out) == 0) {
  cat("\nNo station estimates for either pollutant - nothing written.\n")
  cat("Lower MIN_OBS and MIN_DEATHS at lines 42-43 and re-run.\n")
} else {
  rownames(out) <- NULL
  write.csv(out, "Output/Gap3_regional_heterogeneity.csv", row.names = FALSE)
  cat("Saved: Output/Gap3_regional_heterogeneity.csv\n")
}

#===============================================================================
# RECONCILIATION WITH SECTION 5.3 - READ THIS
#===============================================================================

cat("\n===============================================================\n")
cat("RECONCILIATION WITH SECTION 5.3\n")
cat("===============================================================\n")

pm10 <- regional_all[["PM10AVG"]]
targets <- c("alor setar", "sungai petani", "pegoh ipoh", "batu muda")

if (is.null(pm10) || nrow(pm10) == 0) {
  cat("\nNo PM10 station estimates were produced - skipping the station lookup.\n")
} else {
  cat("\nStations named in Section 5.3 (matched case-insensitively):\n")
  for (t in targets) {
    hit <- pm10[grepl(t, tolower(pm10$Station), fixed = TRUE), ]
    if (nrow(hit) == 0) {
      cat(sprintf("  %-16s NOT FOUND under this name in LOCATION_clean\n", t))
    } else {
      for (i in seq_len(nrow(hit))) {
        cat(sprintf("  %-16s beta = %.5f | %.0f%% OF national | %+.0f%% ABOVE national\n",
                    hit$Station[i], hit$Coefficient[i],
                    hit$Pct_OF_national[i], hit$Pct_ABOVE_national[i]))
      }
    }
  }
}

cat("\n*** ARITHMETIC PROBLEM IN THE CURRENT TEXT ***\n")
cat("Section 5.3 says Kedah stations are '350-380% above the national\n")
cat("average', citing PM10 coefficients of 0.0506-0.0508 against 0.0141.\n")
cat("But 0.0506 / 0.0141 = 3.59, so those stations are approximately 359%\n")
cat("OF the national average, i.e. only about 259% ABOVE it.\n")
cat("'350-380%' is the ratio, not the excess. Either:\n")
cat("   (a) rewrite as 'roughly 3.6 times the national average', or\n")
cat("   (b) keep 'above' and change the figure to ~260%.\n")
cat("The same conflation affects the '42-161% above average' claim for the\n")
cat("urban stations. Both columns are in the CSV so you can pick one\n")
cat("convention and apply it consistently.\n")

cat("\n*** SECOND CAUTION ***\n")
cat("Check the Weak_IV column before quoting any single station. With ~1,270\n")
cat("station-days each, individual first stages are far weaker than the\n")
cat("pooled one, and a station with F<10 should not anchor a paragraph about\n")
cat("a north-south gradient. If several of the named stations are weak, the\n")
cat("honest framing is a correlation between latitude and effect size across\n")
cat("all 45, not four hand-picked stations.\n")

cat("\nGap 3 complete.\n")


