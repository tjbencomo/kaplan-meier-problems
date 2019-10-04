# File: simulations.R
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)

simulate_data <- function(N, maxtime, beta, baseline_hazard, seed) {
  # See Bender et al. 2005 for more information on parametric survival simulation
  # Assume continuous covariate follows a standard normal distribution
  covariate_mean <- 0
  covariate_sd <- 1
  covariate <- rnorm(N, covariate_mean, covariate_sd) # create covariate
  cens <- maxtime * runif(N) # create censoring times (unrelated to survival) --> assumption of Cox
  h <- baseline_hazard * exp(beta * (covariate - covariate_mean)) # create hazard for each patient
  dt <- -log(runif(N)) / h # calculate time to event using exponential with proportional hazards
  e <- ifelse(dt <= cens, 1, 0) # apply censoring
  dt <- pmin(dt, cens) # correct times if censoring applied
  S <- Surv(dt, e)
  linear.fit <- coxph(S ~ covariate)
  median.fit <- coxph(S ~ covariate > median(covariate))
  linear.pvalue <- coef(summary(linear.fit))[, 5]
  median.pvalue <- coef(summary(median.fit))[, 5]
  c(linear = linear.pvalue,
    median = median.pvalue,
    samples = N,
    beta = beta,
    seed = seed)
}

simulate_trial <- function(n_sim, N, maxtime, beta, baseline, seed) {
  set.seed(seed)
  results <- replicate(n_sim, simulate_data(N, maxtime, beta, baseline, seed))
  data.frame(t(results))
}

build_design_matrix <- function(initial_seed, num_seeds, sample_sizes, betas) {
  set.seed(initial_seed)
  seeds <- sample.int(100000, num_seeds)
  design <- expand.grid(
    N = sample_sizes,
    logHR = betas,
    seed = seeds
  )
}

compute_power <- function(simulation_results) {
  power_info <- simulation_results %>%
    gather("linear", "median", 
           key="model", value="pvalue") %>%
    group_by(model, samples, beta, seed) %>%
    summarise(power = sum(pvalue < .05) / n()) %>%
    group_by(model, samples, beta) %>%
    summarise(mean_power = mean(power),
              se = sd(power) / sqrt(n()),
              lower = mean_power - se,
              upper = mean_power + se,
              lower_ci = mean_power - 1.96 * se,
              upper_ci = mean_power + 1.96 * se) %>%
    mutate(HR = paste("HR =", exp(beta)))
}

plot_results <- function(power_results) {
  # Plots 95% confidence intervals, NOT standard error
  power_plot <- power_results %>% 
    ggplot(aes(samples, mean_power)) + 
    geom_point() + 
    # geom_errorbar(aes(ymin = lower, ymax = upper), width=20) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width=20) +
    geom_line(aes(color = model)) + 
    facet_wrap(~HR) +
    geom_hline(yintercept = .8, color="red", linetype="dashed") +
    labs(x = "Samples", y = "Power") +
    guides(color=guide_legend(title="Method"))
}

run_simulation <- function() {
  n_sims <- 10
  maxtime <- 30
  baseline <- .05
  
  sample_sizes <- c(100, 200, 300, 400, 500, 600, 700)
  betas <- log(c(1.0, 1.1, 1.15, 1.2, 1.25, 1.3))
  setup_seed <- 42
  n_seeds <- 10
  
  design <- build_design_matrix(setup_seed, n_seeds, sample_sizes, betas)
  results <- design %>%
    rowwise() %>%
    do(simulate_trial(n_sims, .$N, maxtime, .$logHR, baseline, .$seed))
  power_results <- compute_power(results)
  power_plot <- plot_results(power_results)
  print(power_plot)
}

run_simulation()