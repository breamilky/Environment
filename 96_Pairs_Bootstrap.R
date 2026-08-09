################################################################################
# 96_Pairs_Bootstrap.R
#
# DATE-LEVEL PAIRS CLUSTER BOOTSTRAP  ->  MANUSCRIPT TABLE A2, ROW 1 (PRIMARY)
#
# Target values to reproduce: PM2.5 SE = 0.00368, CI [0.0074, 0.0225]
#                             PM10  SE = 0.00329, CI [0.0068, 0.0202]
#
# WHY THIS SCRIPT EXISTS
#   Table A2 designates the date-level pairs cluster bootstrap as the PRIMARY
#   inference method, and the SE it produces is quoted in the abstract, S4.1,
#   Table 5, Table 11's note and Table A2 itself.
#
#   That estimator is implemented in exactly one place: We_run_first.R line 903,
#   inside run_control_function_iv(). 0202_first.R line 908 implements a
#   DIFFERENT estimator - a moving-block bootstrap over 30/60/90-day blocks -
#   which is what wrote the Bootstrap_SE_30/60/90 columns in your results CSV.
#   Those columns correspond to no row in Table A2.
#
#   This script lifts the pairs bootstrap out of We_run_first.R so it can run
#   on its own, against the same data 0202_first.R used, without re-running
#   that entire script or editing 0202_first.R.
#
# THREE FIXES APPLIED VERSUS THE ORIGINAL
#   1. Validity threshold is now proportional to n_boot. The original hard-coded
#      `n_valid >= 50` in 0202_first.R meant a run with n_boot = 50 needed a
#      100% success rate, which is why blocks 30 and 90 returned NA.
#   2. Truncation tightened from 10x to 3x the point estimate, and BOTH are
#      reported so you can see what the loose rule was admitting.
#   3. The full replicate distribution is saved, so an inflated SE can be
#      traced to specific draws rather than guessed at.
#
# RUNTIME: ~45-90 min per pollutant at n_boot = 500, serial.
#          Set N_CORES > 1 on Mac/Linux to cut this roughly linearly.
#
# Output: Output/TableA2_pairs_bootstrap.csv
#         Output/TableA2_bootstrap_draws.rds
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(sandwich)
})

if (!dir.exists("Output")) dir.create("Output", recursive = TRUE)

#===============================================================================
# CONFIGURATION
#===============================================================================

N_BOOT     <- 10     # Table A2's note says 500; its header says n=1000.
# Pick one and make the manuscript agree - see note at end.
TRUNC_MULT <- 3       # reject draws beyond this multiple of the point estimate
N_CORES    <- 1       # >1 uses parallel::mclapply (Mac/Linux only)
SEED       <- 12345

#===============================================================================
# 1. LOAD DATA
#===============================================================================
# Prefers the RDS written by 0202_first.R so the bootstrap runs on exactly the
# sample the point estimates came from. Falls back to rebuilding from the CSV.

load_pollutant_data <- function(pollutant) {
  
  tag <- if (pollutant == "PM2.5AVG") "pm25" else "pm10"
  
  rds <- c(sprintf("Output/RDS_files/merged_daily_%s_2.rds", tag),
           sprintf("Output/RDS_files/merged_daily_%s.rds",  tag))
  rds <- rds[file.exists(rds)][1]
  
  if (!is.na(rds)) {
    cat("  Loading:", rds, "\n")
    d <- readRDS(rds)
  } else {
    csv <- c("Data/Mid_process_data/merged_death_data_from_script_01.csv",
             "merged_death_data_from_script_01.csv")
    csv <- csv[file.exists(csv)][1]
    if (is.na(csv)) stop("Found neither the RDS nor the merged CSV.")
    cat("  Loading:", csv, "  (RDS not found - rebuilding)\n")
    d <- read.csv(csv, stringsAsFactors = FALSE)
    d <- d %>%
      mutate(date       = as.Date(date),
             station_id = LOCATION_clean,
             year       = lubridate::year(date),
             month      = lubridate::month(date),
             FRP_u_wind_station = total_frp * u_wind,
             FRP_v_wind_station = total_frp * v_wind)
  }
  
  need <- c(pollutant, "deaths", "FRP_u_wind_station", "FRP_v_wind_station",
            "RELATIVEHUMIDITYAVG", "AmbientTemperatureAVG",
            "year", "month", "station_id", "date")
  d <- d[stats::complete.cases(d[, need]), , drop = FALSE]
  
  d$station_id <- factor(d$station_id)
  d$month      <- factor(d$month, levels = sort(unique(d$month)))
  d$year       <- factor(d$year,  levels = sort(unique(d$year)))
  d$date       <- as.Date(d$date)
  
  cat(sprintf("    %s rows, %d stations, %d dates, %s deaths\n",
              format(nrow(d), big.mark = ","),
              nlevels(d$station_id), length(unique(d$date)),
              format(sum(d$deaths), big.mark = ",")))
  d
}

#===============================================================================
# 2. POINT ESTIMATE (identical spec to 0202_first.R)
#===============================================================================

fs_rhs <- paste("FRP_u_wind_station + FRP_v_wind_station +",
                "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                "month + year + station_id")

point_estimate <- function(data, pollutant) {
  
  fs <- lm(as.formula(paste(pollutant, "~", fs_rhs)), data = data)
  data$cf_residuals <- residuals(fs)
  
  ss <- glm(as.formula(paste("deaths ~", pollutant, "+ cf_residuals +",
                             "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                             "month + year + station_id")),
            family = poisson(link = "log"), data = data,
            control = glm.control(maxit = 100, epsilon = 1e-6))
  
  list(coef = as.numeric(coef(ss)[pollutant]), fs = fs, ss = ss)
}

#===============================================================================
# 3. THE PAIRS CLUSTER BOOTSTRAP
#===============================================================================
# Sample n_dates dates WITH REPLACEMENT, retain every station observation for
# each drawn date, re-estimate both stages. This is the estimator Table A2
# describes. It is NOT a block bootstrap - no contiguous windows are drawn.

one_replicate <- function(b, data_by_date, unique_dates, pollutant, coef_point) {
  
  tryCatch({
    
    n_dates  <- length(unique_dates)
    idx      <- sample(seq_len(n_dates), n_dates, replace = TRUE)
    nm       <- as.character(unique_dates[idx])
    
    boot_data <- do.call(rbind, lapply(seq_along(nm), function(i) {
      d <- data_by_date[[nm[i]]]
      d$boot_time_id <- i
      d
    }))
    
    if (nrow(boot_data) < 30)                          return(c(NA, 1))
    if (length(unique(boot_data$station_id)) < 5)       return(c(NA, 2))
    if (length(unique(boot_data$year)) < 2)             return(c(NA, 3))
    
    boot_data$station_id <- factor(boot_data$station_id)
    boot_data$month      <- factor(boot_data$month)
    boot_data$year       <- factor(boot_data$year)
    
    fs_b <- lm(as.formula(paste(pollutant, "~", fs_rhs)), data = boot_data)
    boot_data$cf_resid_boot <- residuals(fs_b)
    
    ss_b <- glm(as.formula(paste("deaths ~", pollutant, "+ cf_resid_boot +",
                                 "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                                 "month + year + station_id")),
                family = poisson(link = "log"), data = boot_data,
                control = glm.control(maxit = 100, epsilon = 1e-6))
    
    if (!pollutant %in% names(coef(ss_b)))              return(c(NA, 4))
    cb <- as.numeric(coef(ss_b)[pollutant])
    if (is.na(cb) || is.infinite(cb))                   return(c(NA, 5))
    
    # Return the raw draw plus a truncation flag. Truncation is applied AFTER
    # the loop so both rules can be compared on the same set of draws.
    c(cb, 0)
    
  }, error = function(e) c(NA, 6))
}

run_pairs_bootstrap <- function(data, pollutant, n_boot = N_BOOT) {
  
  cat("\n  Point estimate...\n")
  pe <- point_estimate(data, pollutant)
  coef_point <- pe$coef
  cat(sprintf("    beta = %.6f\n", coef_point))
  
  cat("  Pre-splitting by date...\n")
  data_by_date <- split(data, data$date)
  unique_dates <- sort(unique(data$date))
  
  cat(sprintf("  Running %d pairs-cluster replicates...\n", n_boot))
  set.seed(SEED)
  t0 <- Sys.time()
  
  if (N_CORES > 1 && .Platform$OS.type != "windows") {
    res <- parallel::mclapply(seq_len(n_boot), one_replicate,
                              data_by_date = data_by_date,
                              unique_dates = unique_dates,
                              pollutant    = pollutant,
                              coef_point   = coef_point,
                              mc.cores     = N_CORES)
  } else {
    res <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      res[[b]] <- one_replicate(b, data_by_date, unique_dates,
                                pollutant, coef_point)
      if (b %% 25 == 0) {
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("    %d/%d  (%.1f min elapsed, ~%.0f min remaining)\n",
                    b, n_boot, el, el / b * (n_boot - b)))
      }
    }
  }
  
  draws  <- vapply(res, function(x) x[1], numeric(1))
  reason <- vapply(res, function(x) x[2], numeric(1))
  
  cat(sprintf("  Done in %.1f minutes\n\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  #--- Failure breakdown -------------------------------------------------------
  labs <- c("ok", "too_few_rows", "too_few_stations", "too_few_years",
            "coef_absent", "coef_na", "error")
  cat("  Replicate outcomes:\n")
  for (k in 0:6) {
    n <- sum(reason == k, na.rm = TRUE)
    if (n > 0) cat(sprintf("    %-18s %d\n", labs[k + 1], n))
  }
  
  #--- Apply both truncation rules to the SAME draws ---------------------------
  finite  <- draws[is.finite(draws)]
  
  summarise_rule <- function(mult) {
    keep <- finite[abs(finite) <= abs(coef_point) * mult]
    n_ok <- length(keep)
    if (n_ok < max(20, 0.8 * n_boot)) {
      return(list(mult = mult, n = n_ok, se = NA, lo = NA, hi = NA,
                  mean = NA, dev = NA,
                  note = sprintf("only %d/%d valid - below 80%% threshold",
                                 n_ok, n_boot)))
    }
    list(mult  = mult,
         n     = n_ok,
         se    = sd(keep),
         lo    = unname(quantile(keep, 0.025)),
         hi    = unname(quantile(keep, 0.975)),
         mean  = mean(keep),
         dev   = 100 * abs(mean(keep) - coef_point) / abs(coef_point),
         note  = sprintf("%d truncated at %gx", length(finite) - n_ok, mult))
  }
  
  strict <- summarise_rule(TRUNC_MULT)
  loose  <- summarise_rule(10)
  
  cat("\n  Truncation rule comparison:\n")
  for (r in list(strict, loose)) {
    cat(sprintf("    %2gx : SE = %-10s  CI [%s, %s]  n = %d  %s\n",
                r$mult,
                ifelse(is.na(r$se), "NA", sprintf("%.6f", r$se)),
                ifelse(is.na(r$lo), "NA", sprintf("%.5f", r$lo)),
                ifelse(is.na(r$hi), "NA", sprintf("%.5f", r$hi)),
                r$n, r$note))
  }
  
  if (!is.na(strict$dev) && strict$dev > 5) {
    cat(sprintf("\n  ! Bootstrap mean deviates from point estimate by %.1f%%.\n",
                strict$dev))
    cat("    Table A2's note claims 2.1%. A larger gap means the draws are\n")
    cat("    skewed - inspect Output/TableA2_bootstrap_draws.rds before using.\n")
  }
  
  list(pollutant = pollutant, coef = coef_point,
       strict = strict, loose = loose, draws = draws, reason = reason)
}

#===============================================================================
# 4. RUN
#===============================================================================

cat("################################################################\n")
cat("# DATE-LEVEL PAIRS CLUSTER BOOTSTRAP (Table A2, PRIMARY)\n")
cat("################################################################\n")

results <- list()

for (poll in c("PM2.5AVG", "PM10AVG")) {
  cat("\n===============================================================\n")
  cat(poll, "\n")
  cat("===============================================================\n")
  d <- load_pollutant_data(poll)
  results[[poll]] <- run_pairs_bootstrap(d, poll)
}

#===============================================================================
# 5. OUTPUT
#===============================================================================

tab <- do.call(rbind, lapply(results, function(r) data.frame(
  Pollutant     = r$pollutant,
  Coefficient   = r$coef,
  Bootstrap_SE  = r$strict$se,
  CI_Lower      = r$strict$lo,
  CI_Upper      = r$strict$hi,
  Bootstrap_Mean= r$strict$mean,
  Mean_Dev_Pct  = r$strict$dev,
  N_Valid       = r$strict$n,
  N_Total       = N_BOOT,
  Trunc_Mult    = TRUNC_MULT,
  SE_at_10x     = r$loose$se,
  stringsAsFactors = FALSE
)))
rownames(tab) <- NULL

cat("\n===============================================================\n")
cat("TABLE A2, ROW 1 - PAIRS CLUSTER BOOTSTRAP\n")
cat("===============================================================\n")
print(tab)

write.csv(tab, "Output/TableA2_pairs_bootstrap.csv", row.names = FALSE)
saveRDS(lapply(results, function(r) list(draws = r$draws, reason = r$reason)),
        "Output/TableA2_bootstrap_draws.rds")

cat("\nSaved: Output/TableA2_pairs_bootstrap.csv\n")
cat("Saved: Output/TableA2_bootstrap_draws.rds\n")

cat("\n---------------------------------------------------------------\n")
cat("RECONCILIATION WITH TABLE A2\n")
cat("---------------------------------------------------------------\n")
cat("Manuscript : PM2.5 SE = 0.00368, CI [0.0074, 0.0225]\n")
cat("             PM10  SE = 0.00329, CI [0.0068, 0.0202]\n\n")
cat("If these reproduce, the pairs bootstrap is confirmed as the source of\n")
cat("the manuscript's primary SEs, and the Bootstrap_SE_30/60/90 columns\n")
cat("from 0202_first.R are a separate block-bootstrap sensitivity check that\n")
cat("should either be reported as such or dropped.\n\n")
cat("If they do NOT reproduce, neither implemented bootstrap produced 0.00368\n")
cat("and the provenance of that number needs establishing before submission.\n\n")
cat("SEPARATELY: Table A2's header says n=1000 but its note says 500\n")
cat("replications. This run used", N_BOOT, "- make the manuscript agree.\n")