# kaplan-meier-problems
Code for *Kaplan Meier Mistakes* article on Medium. `simulations.R` is the script that generates simulation data, measures statistical power, and plots results. `cox_simulations.R` is a first draft script that used a different simulation method.
## Simulation Study Parameters
Survival data is simulated to mimic glioblastoma survival rates.
An exponential distribution is used to model survival times with a
constant baseline hazard of 5%, representing an average survival
of 15 months. Subjects surviving longer than 30 months are censored.
Multiple simulations are generated using multiple seeds to account
for simulation uncertainty. The continuous covariate of interest is
linearly related to survival and covariate values are drawn from a standard
normal distribution. Power is computed using an `alpha=.05`. Equations from 
[Bender et al.](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2059)
are used to simulate parametric survival times. Code was written with 
help from [BBR Section 18.3.3](http://hbiostat.org/doc/bbr.pdf).
## Future Work
### Other Distributions
Bender et al describes survival time generation from Exponential, Weibull,
and Gompertz distributions. The Weibull and Gompertz distribution allows
for the generation of more realistic survival distributions with less
stringent assumptions. It would be interesting to assess power for
specific disease settings using these distributions.
### Nonlinear Covariates
Although many statistical methods assume a linear relationship
between explanatory variable and outcome, this assumption is
not always valid. Splines, fractional polynomials, and other methods
are used to model nonlinear relationships. `simulations.R` contains code
to evaluate the statistical power of restricted cubic splines to identify
linear relationships. These results were excluded from the article due to
their statistical complexity. The spline approach yielded better power than
median Kaplan Meier but was inferior to linear Cox regression. Frank Harrell
has a good explanation of these results at 
[Datamethods](https://discourse.datamethods.org/t/simulating-survival-data-to-evaluate-statistical-power-for-survival-analysis/2486). It would be interesting to evaluate statistical power
with a nonlinear continuous covariate and examine how this affects
Kaplan Meier, and linear/nonlinear Cox regressions.
