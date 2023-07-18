library(DSAIRM)
source("simulate_pkpdbacteria_ode.R")

sim <- simulate_pkpdbacteria_ode(B = 100, I = 1, g = 1, Bmax = 1e+05, 
                            dB = 0.5, kI = 1e-4, r = 1e-4, dI = 2,
                            C0 = 1, dC = 1, C50 = 1, k = 1, Emax = 0,
                            txstart = 10, txinterval = 1, tstart = 0, tfinal = 200, dt = 0.01)

generate_ggplot(list(sim))
