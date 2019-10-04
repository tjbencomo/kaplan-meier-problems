# kaplan-meier-problems
Repository with code to reproduce power simulation figure in the Medium article. `cox_simulations.R` is the script that runs the simulation and records the power for each experiment.
## Simulation Study Parameters
Survival data is simulated to mimic glioblastoma survival rates 
using the `simsurv` package. Data is drawn from a Gompertz distribution
with an approximate mean survival time of 15 months and standard deviation
of 6 months. Scale and shape parameters were calculated using formulas derived
by [Bender et al.](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2059)
The explanatory variable is drawn from a standard normal distribution.  
1000 trials were simulated for each combination of sample size and hazard ratio.
Using an alpha of .05, power was computed by counting the number of trials with p-values less than alpha then dividing by the total number of simulations (1000 trials).
