# File: simulations.R
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(rms)

simulate_data <- function(N, maxtime, beta, baseline_hazard, seed) {
  # See Bender et al. 2005 for more information on parametric survival simulation
  # Assume continuous covariate follows a standard normal distribution and
  # has a linear effect on survival
  covariate_mean <- 0
  covariate_sd <- 1
  covariate <- rnorm(N, covariate_mean, covariate_sd) # create covariate
  cens <- maxtime * runif(N) # create censoring times (unrelated to survival)
  h <- baseline_hazard * exp(beta * (covariate - covariate_mean))
  dt <- -log(runif(N)) / h # Bender formula for times T ~ Exponential(baseline_hazard)
  e <- ifelse(dt <= cens, 1, 0)
  dt <- pmin(dt, cens)
  S <- Surv(dt, e)
  linear.fit <- coxph(S ~ covariate)
  median.fit <- coxph(S ~ covariate > median(covariate))
  # nonlinear.fit <- cph(S ~ rcs(covariate, 3))
  linear.pvalue <- coef(summary(linear.fit))[, 5]
  median.pvalue <- coef(summary(median.fit))[, 5]
  # nonlinear.pvalue <- anova(nonlinear.fit)['covariate', 'P']
  c(linear = linear.pvalue,
    median = median.pvalue,
    # nonlinear = nonlinear.pvalue,
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
    gather("linear", "median", #"nonlinear", 
           key="model", value="pvalue") %>%
    mutate(model = factor(model)) %>%
    mutate(model = fct_recode(model,
                              `Linear Cox` = "linear",
                              `Median KM` = "median"
                              # `Nonlinear Cox` = "nonlinear"
                              )) %>%
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
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width=20) +
    geom_line(aes(color = factor(model))) + 
    facet_wrap(~HR) +
    geom_hline(yintercept = .8, color="red", linetype="dashed") +
    labs(x = "Samples", y = "Power") +
    guides(color=guide_legend(title="Method")) +
    scale_x_continuous("Samples", breaks=seq(0, 600, by=100), limits=c(0, 650)) +
    theme(text = element_text(size=20))
}

run_simulation <- function(save_figure=F) {
  n_sims <- 400
  # n_sims <- 5
  maxtime <- 30
  baseline <- .05
  
  sample_sizes <- c(50, 100, 200, 300, 400, 500, 600)
  betas <- log(c(1.0, 1.1, 1.15, 1.2, 1.25, 1.3))
  # sample_sizes <- c(50, 100, 200, 300, 400, 500, 600)
  # betas <- log(c(1.0, 1.1, 1.2))
  setup_seed <- 58
  n_seeds <- 5
  
  design <- build_design_matrix(setup_seed, n_seeds, sample_sizes, betas)
  results <- design %>%
    rowwise() %>%
    do(simulate_trial(n_sims, .$N, maxtime, .$logHR, baseline, .$seed))
  power_results <- compute_power(results)
  power_plot <- plot_results(power_results)
  print(power_plot)
  if(save_figure) {
    ggsave(file.path(getwd(), "power_plot.png"), width = 12, height = 8)
  }
}

run_simulation(save_figure=T)