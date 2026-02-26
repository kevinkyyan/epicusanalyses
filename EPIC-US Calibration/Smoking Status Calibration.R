# Load libraries
suppressPackageStartupMessages({
  library(epicUS)
  library(dplyr)
  library(tidyr)
  library(future.apply)
  library(knitr)
})

# Setup parallel processing
plan(multisession)

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
# To align the model with National Health Interview Survey (NHIS)
# data each individual is assigned a status: Current, Never (conditional),
# or Former.

# NHIS 2018 and 2023 Prevalence and Trend (AAPC) Targets
targets <- list(
  current2018    = 0.132,
  former2018     = 0.293,
  never2018      = 0.575,
  current2023    = 0.097,
  aapc_post      = -0.019
)

# Weights to prioritize specific targets during optimization
weights <- c(current2018=1, former2018=1, never2018=1, current2023=1, aapc_post=1)

# ------------------------------------------------------------------------------
# Define the Objective Function (RMSE)
# ------------------------------------------------------------------------------
# This function calculates the Root Mean Squared Error (RMSE) between EPIC-US 
# simulated smoking trends and comparing them to the validation targets 
# by modifying intercept parameters.

rmse_objective <- function(params, base_input, targets, wts) {
  epicUS::init_session(settings = settings)
  cal_input <- base_input

  # Map parameters to simulation intercepts
  cal_input$values$smoking$logit_p_current_smoker_0_betas[1, "Intercept"] <- params[1]
  cal_input$values$smoking$logit_p_never_smoker_con_not_current_0_betas[1, "intercept"] <- params[2]
  cal_input$values$smoking$ln_h_ces_betas["intercept"] <- params[3]

  # Run EPIC-US
  epicUS::run(input = cal_input$values)
  output_data <- epicUS::Cget_output_ex()
  epicUS::terminate_session()

  # Convert raw counts to prevalence percentages
  status_counts <- output_data$n_smoking_status_by_ctime
  prevalence_df <- as.data.frame(status_counts / rowSums(status_counts))
  colnames(prevalence_df) <- c("Never", "Current", "Former")
  prevalence_df$Year <- 2015:(2015 + nrow(prevalence_df) - 1)

  # Calculate trend (AAPC) for the post-2023 projection period
  future_years <- prevalence_df[prevalence_df$Year %in% 2024:2060, ]
  aapc_simulated <- exp(coef(lm(log(Current) ~ Year, data = future_years))[["Year"]]) - 1

  # Gather simulated trends for comparison
  simulated_metrics <- c(
    current2018    = prevalence_df$Current[prevalence_df$Year == 2018],
    former2018     = prevalence_df$Former[prevalence_df$Year == 2018],
    never2018      = prevalence_df$Never[prevalence_df$Year == 2018],
    current2023    = prevalence_df$Current[prevalence_df$Year == 2023],
    aapc_post      = aapc_simulated
  )

  # Weighted RMSE calculation
  deviations <- unlist(simulated_metrics[names(targets)]) - unlist(targets)
  total_rmse <- sqrt(mean((deviations^2) * unlist(wts[names(deviations)])))

  return(total_rmse)
}

# ------------------------------------------------------------------------------
# Smoking Status Optimization
# ------------------------------------------------------------------------------
# Use nlminb constrained to ±1 around starting values.

initial_guess <- c(b0_cur = -0.2, b0_never = 3.7, b0_ces = -3.7)
lower_bounds  <- initial_guess - 1
upper_bounds  <- initial_guess + 1

message("General search for different combination of parameters")
coarse_grid <- expand.grid(
  b0_cur = seq(lower_bounds[1], upper_bounds[1], length.out = 5),
  b0_never = seq(lower_bounds[2], upper_bounds[2], length.out = 5),
  b0_ces = seq(lower_bounds[3], upper_bounds[3], length.out = 5)
)

grid_results <- future_apply(coarse_grid, 1, function(row) {
  rmse_objective(as.numeric(row), base_input=input, targets=targets, wts=weights)
}, future.seed = TRUE)

refined_start <- as.numeric(coarse_grid[which.min(grid_results), 1:3])

message("Finding parameters using nlminb")
result_nlminb <- nlminb(
  start = refined_start,
  objective = rmse_objective,
  base_input = input,
  targets = targets,
  wts = weights,

  # Bounds are still required HERE for the optimizer
  lower = lower_bounds,
  upper = upper_bounds,

  # CLEANED: Removed 'lower_bounds' and 'upper_bounds' passed to function

  control = list(eval.max = 100, iter.max = 100)
)

final_params <- result_nlminb$par

# ------------------------------------------------------------------------------
# Final Adjustments & Verification
# ------------------------------------------------------------------------------

input$values$smoking$logit_p_current_smoker_0_betas[1, "Intercept"] <- final_params[1]
input$values$smoking$logit_p_never_smoker_con_not_current_0_betas[1, "intercept"] <- final_params[2]
input$values$smoking$ln_h_ces_betas["intercept"] <- final_params[3]
input$values$smoking$pack_years_0_betas[1, "intercept"] <- 30

epicUS::init_session(settings = settings)
epicUS::run(input = input$values)
verification_output <- epicUS::Cget_output_ex()
epicUS::terminate_session()

# ------------------------------------------------------------------------------
# Diagnostics Table
# ------------------------------------------------------------------------------
final_counts <- verification_output$n_smoking_status_by_ctime
final_prevalence <- as.data.frame(final_counts / rowSums(final_counts))
colnames(final_prevalence) <- c("Never","Current","Former")
final_prevalence$Year <- 2015:(2015 + nrow(final_prevalence) - 1)

diag_tab <- data.frame(
  Metric = c("Current 2018", "Former 2018", "Never 2018", "Current 2023"),
  Target = c(targets$current2018, targets$former2018, targets$never2018, targets$current2023),
  Simulated = c(final_prevalence$Current[final_prevalence$Year==2018],
                final_prevalence$Former[final_prevalence$Year==2018],
                final_prevalence$Never[final_prevalence$Year==2018],
                final_prevalence$Current[final_prevalence$Year==2023])
)

print(kable(diag_tab, digits = 4, caption = "Final Smoking Calibration Results"))
