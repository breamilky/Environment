################################################################################
# 90_Common_Setup.R
#
# SHARED SETUP FOR THE FIVE "MISSING ANALYSIS" SCRIPTS (91-95)
#
# Purpose : Load data, build station-specific instruments, and define ONE
#           canonical CF-IV estimator so all five scripts produce numbers on
#           exactly the same basis as the manuscript's main specification.
#
# Spec replicated here (manuscript Sections 3.2-3.3):
#   First stage : Pollutant ~ FRP_u_wind_station + FRP_v_wind_station
#                             + RH + Temp + month FE + year FE + station FE
#   Second stage: deaths ~ OBSERVED Pollutant + first-stage residual
#                          + RH + Temp + month FE + year FE + station FE
#                 estimated by Poisson (or Negative Binomial in script 92)
#   Inference   : two-way clustered by station and date (Cameron, Gelbach
#                 and Miller 2011)
#
# DO NOT RUN THIS FILE ON ITS OWN. Each of 91-95 sources it at the top.
#
# Run order (nothing here depends on the 01_Feb2 / 0202 / We_run RDS files,
# so these five can be run in any order, in parallel, straight after
# 01_Clean_Merge.R):
#
#   source("91_Exclude2020.R")            -> Table 11, COVID row
#   source("92_NegativeBinomial.R")       -> Table 11, NegBin row
#   source("93_Regional_Heterogeneity.R") -> Section 5.3
#   source("94_FireDays_Above5.R")        -> Table 12, final column
#   source("95_Welfare_Sensitivity.R")    -> Section 6.7
################################################################################


#===============================================================================
# 0. PACKAGES
#===============================================================================

.pkgs <- c("dplyr", "tidyr", "lubridate", "sandwich", "lmtest", "MASS")

for (p in .pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(sandwich)
  library(lmtest)
})
# NOTE: MASS is called as MASS::glm.nb rather than attached, because
# library(MASS) masks dplyr::select and will silently break the pipelines below.

if (!dir.exists("Output")) dir.create("Output", recursive = TRUE)


#===============================================================================
# 1. LOAD AND PREPARE DATA
#===============================================================================

prep_analysis_data <- function(
    paths = c("Data/Mid_process_data/merged_death_data_from_script_01.csv",
              "merged_death_data_from_script_01.csv")) {
  
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) {
    stop("Could not find merged_death_data_from_script_01.csv. ",
         "Run 01_Clean_Merge.R first, or edit the `paths` argument.")
  }
  cat("Loading:", path, "\n")
  
  d <- read.csv(path, stringsAsFactors = FALSE)
  
  d <- d %>%
    mutate(
      date       = as.Date(date),
      station_id = as.factor(LOCATION_clean),
      year       = lubridate::year(date),
      month      = lubridate::month(date)
    )
  
  # --- Station-specific wind vectors and instruments -------------------------
  # u = speed * sin(direction), v = speed * cos(direction)  [manuscript S2.2]
  # Rebuilt here rather than trusted from the CSV, because the older scripts
  # (02-06) wrote a NON-station-specific FRP_u_wind / FRP_v_wind under similar
  # names. Overwriting guarantees the correct instrument regardless of which
  # upstream script produced the file.
  d <- d %>%
    mutate(
      FRP_u_wind_station = total_frp * u_wind,
      FRP_v_wind_station = total_frp * v_wind
    )
  
  cat("  Observations :", format(nrow(d), big.mark = ","), "\n")
  cat("  Stations     :", length(unique(d$station_id)), "\n")
  cat("  Dates        :", length(unique(d$date)), "\n")
  cat("  Date range   :", as.character(min(d$date, na.rm = TRUE)), "to",
      as.character(max(d$date, na.rm = TRUE)), "\n")
  cat("  Total deaths :", format(sum(d$deaths, na.rm = TRUE), big.mark = ","), "\n\n")
  
  d
}


#' Restrict to complete cases for a given pollutant.
#' Guarantees nobs(model) == nrow(data), which the clustering code relies on.
complete_sample <- function(data, pollutant) {
  need <- c(pollutant, "deaths", "FRP_u_wind_station", "FRP_v_wind_station",
            "RELATIVEHUMIDITYAVG", "AmbientTemperatureAVG",
            "year", "month", "station_id", "date")
  data[stats::complete.cases(data[, need]), , drop = FALSE]
}


#===============================================================================
# 2. TWO-WAY CLUSTERED VARIANCE
#===============================================================================

#' Two-way clustered VCOV: V_station + V_date - V_intersection.
#' Falls back to station-only if the result is not positive semi-definite,
#' which is the same guard used in We_run_first.R.
twoway_vcov <- function(model, cl_station, cl_date, quiet = FALSE) {
  
  v_s  <- sandwich::vcovCL(model, cluster = cl_station, type = "HC1")
  v_d  <- sandwich::vcovCL(model, cluster = cl_date,    type = "HC1")
  v_sd <- sandwich::vcovCL(model,
                           cluster = paste(cl_station, cl_date, sep = "_"),
                           type = "HC1")
  
  v <- v_s + v_d - v_sd
  
  # Cameron, Gelbach & Miller (2011) eigenvalue correction. The raw sum is
  # routinely not PSD here because V_station has rank <= 44 (45 stations)
  # against ~62 parameters. Zeroing negative eigenvalues gives the nearest
  # PSD matrix; falling back to V_station discards the date dimension
  # entirely, which is why every SE in this run was station-clustered.
  e <- eigen(v, symmetric = TRUE)
  if (any(e$values < -1e-10)) {
    if (!quiet) {
      cat(sprintf("    two-way VCOV: %d negative eigenvalue(s) zeroed (CGM 2011)\n",
                  sum(e$values < -1e-10)))
    }
    v <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
    dimnames(v) <- dimnames(v_s)
  }
  v
}


#===============================================================================
# 3. CANONICAL CF-IV ESTIMATOR
#===============================================================================

#' Control-function IV estimator (manuscript equations 1 and 2).
#'
#' @param data       analysis frame, already complete-cased
#' @param pollutant  "PM2.5AVG" or "PM10AVG"
#' @param family     "poisson" (default) or "negbin"
#' @param use_station_fe  FALSE for single-station regressions (script 93)
#' @param cluster    "twoway" (default), "station", or "none"
#' @param label      printed tag
#'
#' @return list with coefficient, se, p-value, F-statistics, N, and both models
run_cf_iv <- function(data,
                      pollutant,
                      family          = c("poisson", "negbin"),
                      use_station_fe  = TRUE,
                      cluster         = c("twoway", "station", "none"),
                      label           = "",
                      quiet           = FALSE) {
  
  family  <- match.arg(family)
  cluster <- match.arg(cluster)
  
  data <- complete_sample(data, pollutant)
  if (nrow(data) < 50) {
    return(list(ok = FALSE, reason = "fewer than 50 usable observations",
                label = label, n = nrow(data)))
  }
  
  #--- Build FE terms, dropping any that have no variation -------------------
  fe <- character(0)
  if (length(unique(data$month))      > 1) fe <- c(fe, "factor(month)")
  if (length(unique(data$year))       > 1) fe <- c(fe, "factor(year)")
  if (use_station_fe &&
      length(unique(data$station_id)) > 1) fe <- c(fe, "factor(station_id)")
  
  ctrl  <- c("RELATIVEHUMIDITYAVG", "AmbientTemperatureAVG")
  inst  <- c("FRP_u_wind_station", "FRP_v_wind_station")
  rhs   <- paste(c(ctrl, fe), collapse = " + ")
  
  #--- STAGE 1 ---------------------------------------------------------------
  f_full <- as.formula(paste(pollutant, "~", paste(inst, collapse = " + "), "+", rhs))
  f_rest <- as.formula(paste(pollutant, "~", rhs))
  
  fs_full <- try(lm(f_full, data = data), silent = TRUE)
  fs_rest <- try(lm(f_rest, data = data), silent = TRUE)
  if (inherits(fs_full, "try-error") || inherits(fs_rest, "try-error")) {
    return(list(ok = FALSE, reason = "first stage failed", label = label))
  }
  
  data$cf_resid <- residuals(fs_full)
  
  # Partial F for the two excluded instruments (conventional / unclustered).
  # This is the F reported in manuscript Table 3 "unclustered" and Table 11.
  n_obs   <- nobs(fs_full)
  k_par   <- length(coef(fs_full))
  rss_f   <- sum(residuals(fs_full)^2)
  rss_r   <- sum(residuals(fs_rest)^2)
  f_uncl  <- ((rss_r - rss_f) / 2) / (rss_f / (n_obs - k_par))
  
  # Cluster-robust Wald F for the same two instruments (Table 3 "clustered")
  f_clust <- NA_real_
  if (cluster != "none") {
    v_fs <- if (cluster == "twoway") {
      twoway_vcov(fs_full, as.character(data$station_id),
                  as.character(data$date), quiet = TRUE)
    } else {
      sandwich::vcovCL(fs_full, cluster = as.character(data$station_id), type = "HC1")
    }
    b_i <- coef(fs_full)[inst]
    v_i <- v_fs[inst, inst, drop = FALSE]
    f_clust <- try(as.numeric(t(b_i) %*% solve(v_i) %*% b_i) / 2, silent = TRUE)
    if (inherits(f_clust, "try-error")) f_clust <- NA_real_
  }
  
  r2_full   <- summary(fs_full)$r.squared
  r2_rest   <- summary(fs_rest)$r.squared
  partial_r2 <- (r2_full - r2_rest) / (1 - r2_rest)
  
  #--- STAGE 2 ---------------------------------------------------------------
  # NOTE: OBSERVED pollutant + residual, not the fitted value.
  f_ss <- as.formula(paste("deaths ~", pollutant, "+ cf_resid +", rhs))
  
  ss <- try({
    if (family == "poisson") {
      glm(f_ss, family = poisson(link = "log"), data = data)
    } else {
      MASS::glm.nb(f_ss, data = data, control = glm.control(maxit = 100))
    }
  }, silent = TRUE)
  
  if (inherits(ss, "try-error")) {
    return(list(ok = FALSE, reason = paste("second stage failed:", ss),
                label = label))
  }
  
  #--- Inference -------------------------------------------------------------
  vc <- switch(
    cluster,
    twoway  = twoway_vcov(ss, as.character(data$station_id),
                          as.character(data$date), quiet = quiet),
    station = sandwich::vcovCL(ss, cluster = as.character(data$station_id),
                               type = "HC1"),
    none    = vcov(ss)
  )
  
  b  <- coef(ss)[pollutant]
  se <- sqrt(vc[pollutant, pollutant])
  z  <- b / se
  p  <- 2 * pnorm(-abs(z))
  
  # Endogeneity test: is the control-function term significant?
  se_cf <- sqrt(vc["cf_resid", "cf_resid"])
  p_cf  <- 2 * pnorm(-abs(coef(ss)["cf_resid"] / se_cf))
  
  out <- list(
    ok            = TRUE,
    label         = label,
    pollutant     = pollutant,
    family        = family,
    n             = nrow(data),
    n_stations    = length(unique(data$station_id)),
    n_dates       = length(unique(data$date)),
    deaths        = sum(data$deaths),
    coefficient   = unname(b),
    se            = unname(se),
    z             = unname(z),
    p_value       = unname(p),
    ci_lo         = unname(b - 1.96 * se),
    ci_hi         = unname(b + 1.96 * se),
    effect_pct10  = (exp(unname(b) * 10) - 1) * 100,
    ci_lo_pct10   = (exp(unname(b - 1.96 * se) * 10) - 1) * 100,
    ci_hi_pct10   = (exp(unname(b + 1.96 * se) * 10) - 1) * 100,
    f_unclustered = unname(f_uncl),
    f_clustered   = unname(f_clust),
    partial_r2    = unname(partial_r2),
    coef_u        = unname(coef(fs_full)["FRP_u_wind_station"]),
    coef_v        = unname(coef(fs_full)["FRP_v_wind_station"]),
    p_endogeneity = unname(p_cf),
    theta         = if (family == "negbin") ss$theta else NA_real_,
    theta_se      = if (family == "negbin") ss$SE.theta else NA_real_,
    first_stage   = fs_full,
    second_stage  = ss
  )
  
  if (!quiet) print_cf(out)
  out
}


print_cf <- function(r) {
  if (!isTRUE(r$ok)) {
    cat(sprintf("  [%s] FAILED: %s\n", r$label, r$reason)); return(invisible())
  }
  cat(sprintf("  %-34s beta = %9.6f  SE = %8.6f  p = %.4f\n",
              r$label, r$coefficient, r$se, r$p_value))
  cat(sprintf("  %-34s effect/10 ug/m3 = %6.2f%%  [%.2f%%, %.2f%%]\n",
              "", r$effect_pct10, r$ci_lo_pct10, r$ci_hi_pct10))
  cat(sprintf("  %-34s first-stage F: %.0f (uncl.) / %.1f (clust.)   N = %s\n\n",
              "", r$f_unclustered, r$f_clustered, format(r$n, big.mark = ",")))
}


#' Flatten a run_cf_iv result into one row for rbind-ing into a table.
cf_row <- function(r) {
  if (!isTRUE(r$ok)) {
    return(data.frame(Label = r$label, Coefficient = NA, SE = NA, P_value = NA,
                      Effect_per10 = NA, CI_lo_pct = NA, CI_hi_pct = NA,
                      F_unclustered = NA, F_clustered = NA, N = NA,
                      Status = paste("FAILED:", r$reason),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    Label         = r$label,
    Coefficient   = r$coefficient,
    SE            = r$se,
    P_value       = r$p_value,
    Effect_per10  = r$effect_pct10,
    CI_lo_pct     = r$ci_lo_pct10,
    CI_hi_pct     = r$ci_hi_pct10,
    F_unclustered = r$f_unclustered,
    F_clustered   = r$f_clustered,
    N             = r$n,
    Status        = "OK",
    stringsAsFactors = FALSE
  )
}

cat("90_Common_Setup.R loaded: prep_analysis_data(), run_cf_iv(), twoway_vcov()\n\n")

