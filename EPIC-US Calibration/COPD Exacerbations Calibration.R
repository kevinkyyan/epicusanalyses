# Load libraries
suppressPackageStartupMessages({
  library(epicUS)
  library(dplyr)
  library(tidyr)
  library(knitr)
  library(future.apply)
})

# Setup parallel processing
plan(multisession)

# EPIC-US model initialization
settings <- epicUS::get_default_settings()
settings$record_mode   <- 2
settings$n_base_agents <- 3.5e6
epicUS::init_session(settings = settings)
input <- epicUS::get_input()
time_horizon <- 46
input$values$global_parameters$time_horizon <- time_horizon
epicUS::terminate_session()

# Number of EPIC-US runs
N_REPLICATES <- 100

# ------------------------------------------------------------------------------
# Define Calibration Targets
# ------------------------------------------------------------------------------
# To align the model with published  align the model’s outputs with validation 
# targets for COPD exacerbations from Wallace et al. 2019 and 
# Hoogendorn et al. 2021

targets <- list(
  frequency = c(gold1 = 0.82, gold2 = 1.17, gold3 = 1.61, gold4 = 2.10), # Hoogendoorn 2021
  ford_hosp = 510,                                                     # Ford 2015
  mod_sev   = c(gold1 = 0.404, gold2 = 0.489, gold3 = 0.836, gold4 = 0.891), # Wallace 2019
  sev_only  = c(gold1 = 0.120, gold2 = 0.139, gold3 = 0.254, gold4 = 0.422)  # Wallace 2019
)

# ------------------------------------------------------------------------------
# Calculation Function (Parallelized)
# ------------------------------------------------------------------------------
# Processes pooled exacerbation rates across N_REPLICATES using parallel processing
get_pooled_metrics <- function(input) {

  replicate_results <- future_lapply(1:N_REPLICATES, function(i) {

    # Initialize and run session inside the worker
    epicUS::init_session(settings = settings)
    epicUS::run(input = input$values)
    events <- as.data.frame(epicUS::Cget_all_events_matrix())
    output_ex <- epicUS::Cget_output_ex()
    epicUS::terminate_session()

    # Hospitalization rate calculation
    pop_alive <- rowSums(output_ex$n_alive_by_ctime_sex)
    # Rate per 100,000 in 2016
    hosp_i <- (output_ex$n_exac_by_ctime_severity[2, 3] + output_ex$n_exac_by_ctime_severity[2, 4]) * (100000 / pop_alive[2])

    # Initialize counters for this single run
    res <- list(pt = rep(0, 4), tot = rep(0, 4), ms = rep(0, 4), sev = rep(0, 4), hosp = hosp_i)

    # Process event-level data
    dx <- subset(events, diagnosis > 0 & gold > 0)
    if (nrow(dx) > 0) {
      dx <- dx[order(dx$id, dx$local_time), ]

      # Calculate duration
      dur <- pmax(lead(dx$local_time, default = dx$local_time[nrow(dx)]) - dx$local_time, 0)
      dur[dx$id != lead(dx$id, default = dx$id[nrow(dx)])] <- 0

      # Aggregate counts by GOLD stage
      res$pt  <- as.numeric(tapply(dur, factor(dx$gold, 1:4), sum, na.rm=TRUE))
      res$tot <- as.numeric(table(factor(subset(dx, event == 5)$gold, 1:4)))
      res$ms  <- as.numeric(table(factor(subset(dx, event == 5 & exac_status %in% 2:4)$gold, 1:4)))
      res$sev <- as.numeric(table(factor(subset(dx, event == 5 & exac_status %in% 3:4)$gold, 1:4)))

      res$pt[is.na(res$pt)] <- 0
    }
    return(res)
  }, future.seed = TRUE)

  # Pool results across all replicates
  pooled <- Reduce(function(a, b) {
    list(
      pt   = a$pt + b$pt,
      tot  = a$tot + b$tot,
      ms   = a$ms + b$ms,
      sev  = a$sev + b$sev,
      hosp = a$hosp + b$hosp
    )
  }, replicate_results)

  # Calculate Final Rates
  list(
    hosp_2016 = pooled$hosp / N_REPLICATES,
    freq      = pooled$tot / pmax(pooled$pt, 1e-9),
    mod_sev   = pooled$ms  / pmax(pooled$pt, 1e-9),
    sev_only  = pooled$sev / pmax(pooled$pt, 1e-9)
  )
}

# ------------------------------------------------------------------------------
# Total Exacerbation Optimization
# ------------------------------------------------------------------------------
# This function calculates the Root Mean Squared Error (RMSE) between EPIC-US 
# simulated smoking trends and comparing them to the validation targets by 
# modifying intercept parameters and GOLD stage coefficients 

calculate_rmse_frequency <- function(params, input, lower_bounds, upper_bounds) {
  params <- pmax(pmin(params, upper_bounds), lower_bounds)
  input$values$exacerbation$ln_rate_betas[1, c(1, 6:9)] <- params

  sim <- get_pooled_metrics(input)
  rmse <- sqrt(mean((sim$freq - targets$frequency)^2))

  message(sprintf("Freq RMSE: %.4f | Intercept: %.2f", rmse, params[1]))
  return(rmse)
}

# Use nlminb constrained to ±1 around starting values.
initial_guess <- c(intercept = 1.4, g1 = 0.3, g2 = -0.3, g3 = 0.08, g4 = -0.35)
lower_bounds  <- initial_guess - 1
upper_bounds  <- initial_guess + 1

fit_frequency <- nlminb(
  start        = initial_guess,
  objective    = calculate_rmse_frequency,
  input        = input,
  lower        = lower_bounds,
  upper        = upper_bounds,
  lower_bounds = lower_bounds,
  upper_bounds = upper_bounds,
  control      = list(eval.max = 250, iter.max = 250)
)

# ------------------------------------------------------------------------------
# Severity of Exacerbation Optimization
# ------------------------------------------------------------------------------
# This function calculates the Root Mean Squared Error (RMSE) between EPIC-US 
# simulated smoking trends and comparing them to the validation targets by 
# modifying intercept parameters

calculate_rmse_severity <- function(params, input, lower_bounds, upper_bounds, best_freq) {
  params <- pmax(pmin(params, upper_bounds), lower_bounds)

  input$values$exacerbation$ln_rate_betas[1, c(1, 6:9)]  <- best_freq
  input$values$exacerbation$logit_severity_betas[1, 1:3] <- params

  sim <- get_pooled_metrics(input)

  rmse_ms   <- sqrt(mean((sim$mod_sev - targets$mod_sev)^2))
  rmse_sev  <- sqrt(mean((sim$sev_only - targets$sev_only)^2))
  penalty   <- (log(sim$hosp_2016 / targets$ford_hosp))^2

  return(rmse_ms + rmse_sev + penalty)
}

# Use nlminb constrained to ±5 around starting values.
initial_guess_sev <- c(i1 = -3.609, i2 = 2.202, i3 = 5.208)
lower_bounds_sev  <- initial_guess_sev - 5
upper_bounds_sev  <- initial_guess_sev + 5

fit_severity <- nlminb(
  start        = initial_guess_sev,
  objective    = calculate_rmse_severity,
  input        = input,
  best_freq    = fit_frequency$par,
  lower        = lower_bounds_sev,
  upper        = upper_bounds_sev,
  lower_bounds = lower_bounds_sev,
  upper_bounds = upper_bounds_sev,
  control      = list(eval.max = 300, iter.max = 300)
)

# ------------------------------------------------------------------------------
# Final Verification & Results
# ------------------------------------------------------------------------------
input$values$exacerbation$ln_rate_betas[1, c(1, 6:9)]  <- fit_frequency$par
input$values$exacerbation$logit_severity_betas[1, 1:3] <- fit_severity$par

final <- get_pooled_metrics(input)

cat("Final Exacerbation Calibration Results")
kable(data.frame(
  Metric    = c("Ford Hosp 2016", "Total Rate", "Severity"),
  Target    = c(targets$ford_hosp, targets$frequency[1], targets$sev_only[4]),
  Simulated = c(final$hosp_2016, final$freq[1], final$sev_only[4])
), digits = 3)
