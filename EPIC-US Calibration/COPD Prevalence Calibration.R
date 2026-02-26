# Load libraries
suppressPackageStartupMessages({
  library(epicUS)
  library(dplyr)
  library(tidyr)
  library(knitr)
})

# EPIC-US model initialization
settings <- epicUS::get_default_settings()
settings$record_mode   <- 0
settings$n_base_agents <- 1e6
epicUS::init_session(settings = settings)

input <- epicUS::get_input()
time_horizon <- 46
input$values$global_parameters$time_horizon <- time_horizon

# ------------------------------------------------------------------------------
# Define Calibration Targets
# ------------------------------------------------------------------------------
# We recalibrate the sex-specific intercept terms to match NHANES age-specific
# prevalence (40–59 and 60–79 years) while maintaining the 1.4:1 male-to-female
# COPD prevalence ratio observed in NHANES.

# Calibration targets and ratio constraints based on NHANES data
age_targets  <- c(p40_59 = 0.081, p60_79 = 0.144)
ratio_target <- 1.4

# ------------------------------------------------------------------------------
# Define the Objective Function (RMSE)
# ------------------------------------------------------------------------------
# This function calculates the Root Mean Squared Error (RMSE) between EPIC-US 
# simulated prevalence estimates and comparing them to the validation targets by 
# modifying sex specific intercept parameters. 

rmse_objective <- function(params, base_input, targets, ratio_goal, lower_bounds, upper_bounds) {

  # Ensure the optimizer stays within specified bounds
  params <- pmax(pmin(params, upper_bounds), lower_bounds)

  epicUS::init_session(settings = settings)
  cal_input <- base_input

  # Map parameters to COPD intercepts: male and female
  copd_intercepts <- cal_input$values$COPD$logit_p_COPD_betas_by_sex
  copd_intercepts["intercept", "male"]   <- params[1]
  copd_intercepts["intercept", "female"] <- params[2]
  cal_input$values$COPD$logit_p_COPD_betas_by_sex <- copd_intercepts

  # Run EPIC-US
  epicUS::run(input = cal_input$values)
  output_data <- epicUS::Cget_output_ex()
  epicUS::terminate_session()

  # Index for the first year of simulation (Baseline)
  baseline_year <- 1

  # Calculate Age-group prevalence for the baseline year
  prev_40_59 <- sum(output_data$n_COPD_by_ctime_age[baseline_year, 40:59]) /
    sum(output_data$n_alive_by_ctime_age[baseline_year, 40:59])

  prev_60_79 <- sum(output_data$n_COPD_by_ctime_age[baseline_year, 60:79]) /
    sum(output_data$n_alive_by_ctime_age[baseline_year, 60:79])

  # Calculate Sex prevalence for ratio check
  male_prev   <- output_data$n_COPD_by_ctime_sex[baseline_year, 1] /
    output_data$n_alive_by_ctime_sex[baseline_year, 1]
  female_prev <- output_data$n_COPD_by_ctime_sex[baseline_year, 2] /
    output_data$n_alive_by_ctime_sex[baseline_year, 2]
  sim_ratio   <- as.numeric(male_prev / female_prev)

  # Calculate RMSE for age targets
  rmse_age      <- sqrt(mean((c(prev_40_59, prev_60_79) - targets)^2))
  ratio_penalty <- (log(sim_ratio / ratio_goal))^2

  return(rmse_age + ratio_penalty)
}

# ------------------------------------------------------------------------------
# COPD Prevalence Optimization
# ------------------------------------------------------------------------------
# Use nlminb constrained to ±1 around starting values.

initial_guess <- c(male_intercept = -4.522189, female_intercept = -4.074861)
lower_bounds  <- initial_guess - 1
upper_bounds  <- initial_guess + 1

message("General search for different combination of parameter")
grid_search <- expand.grid(m = seq(lower_bounds[1], upper_bounds[1], length.out = 3),
                           f = seq(lower_bounds[2], upper_bounds[2], length.out = 3))

grid_search$score <- apply(grid_search, 1, function(row) {
  rmse_objective(as.numeric(row), base_input = input, targets = age_targets,
                 ratio_goal = ratio_target, lower_bounds = lower_bounds,
                 upper_bounds = upper_bounds)
})
refined_start <- as.numeric(grid_search[which.min(grid_search$score), 1:2])

message("Finding parameters using nlminb")
# Using nlminb for a bounded, gradient-based nonlinear optimization ±1 
# around starting values.
optimization_result <- nlminb(
  start = refined_start,
  objective = rmse_objective,
  base_input = input,
  targets = age_targets,
  ratio_goal = ratio_target,
  lower_bounds = lower_bounds,
  upper_bounds = upper_bounds,
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(eval.max = 500, iter.max = 500)
)

final_params <- optimization_result$par

# ------------------------------------------------------------------------------
# Final Verification & Diagnostics
# ------------------------------------------------------------------------------
input$values$COPD$logit_p_COPD_betas_by_sex["intercept", "male"]   <- final_params[1]
input$values$COPD$logit_p_COPD_betas_by_sex["intercept", "female"] <- final_params[2]

epicUS::init_session(settings = settings)
epicUS::run(input = input$values)
verification_output <- epicUS::Cget_output_ex()
epicUS::terminate_session()

# Table
baseline_year <- 1
diagnostics_table <- data.frame(
  Metric    = c("Age 40–59 (LLN)", "Age 60–79 (LLN)", "Male:Female ratio"),
  Target    = c(age_targets["p40_59"], age_targets["p60_79"], ratio_target),
  Simulated = c(
    sum(verification_output$n_COPD_by_ctime_age[baseline_year, 40:59]) /
      sum(verification_output$n_alive_by_ctime_age[baseline_year, 40:59]),
    sum(verification_output$n_COPD_by_ctime_age[baseline_year, 60:79]) /
      sum(verification_output$n_alive_by_ctime_age[baseline_year, 60:79]),
    (verification_output$n_COPD_by_ctime_sex[baseline_year, 1] / verification_output$n_alive_by_ctime_sex[baseline_year, 1]) /
      (verification_output$n_COPD_by_ctime_sex[baseline_year, 2] / verification_output$n_alive_by_ctime_sex[baseline_year, 2])
  )
)

print(kable(diagnostics_table, digits = 4, caption = "Final COPD Prevalence Calibration Results"))
