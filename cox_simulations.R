library(tidyr)
library(ggplot2)
library(survminer)
library(survival)
library(simsurv)
library(dplyr)
library(readr)

simulate_data <- function(N, maxtime, shape, scale, log_hazard_ratio) {
      covariates <- data.frame(id = 1:N,
                               x = rnorm(N, 0, 1)
                                                                                     )
  times <- simsurv(dist = "gompertz",
                   lambdas = shape,
                   gammas = scale,
                   betas = c(x = log_hazard_ratio),
                   x = covariates,
                   maxt = maxtime)

    df <- merge(covariates, times)
    continuous.fit <- coxph(Surv(eventtime, status) ~ x, data=df)
    median.fit <- coxph(Surv(eventtime, status) ~ x > median(x), data=df)
    continuous.pvalue <- coef(summary(continuous.fit))[, 5]
    median.pvalue <- coef(summary(median.fit))[, 5]
    c(linear = continuous.pvalue,
      median = median.pvalue,
      samples = N,
      beta = log_hazard_ratio)
}

simulate_trial <- function(num_simulations, N, maxtime, shape, scale, effect) {
      results <- replicate(num_sims, simulate_data(N, maxtime, 
                                                   shape, scale, 
                                                   effect))
  data.frame(t(results))
}

sample_sizes <- c(50, 100, 200, 300, 400, 500, 1000)
effect_sizes <- log(c(1.0, 1.05, 1.1, 1.2, 1.3, 1.5))
design <- expand.grid(N = sample_sizes,
                      logHR = effect_sizes)
num_sims <- 1000
mean_survival <- 15 #in months
sd_survival <- 5
maxtime <- 30
# Parameter approximations from Bender et al 
# "Generating Survival Times to Simulate Cox Proportional Hazards Models"
shape <- pi / (sqrt(6) * sd_survival)
scale <- shape * exp(-0.5772 - (shape * mean_survival))

set.seed(56)
results <- design %>%
  rowwise() %>%
  do(simulate_trial(num_simulations, .$N, maxtime, scale, shape, .$logHR))

power_info <- results %>%
  gather("linear", "median", 
                  key="model", value="pvalue") %>%
  group_by(model, samples, beta,) %>%
    summarise(power = sum(pvalue < .05) / n()) %>%
    mutate(HR = paste("HR =", exp(beta)))
save_dir <- '/scratch/users/tbencomo'
write_csv(power_info, file.path(save_dir, 'cox_power_simulations.csv'))
