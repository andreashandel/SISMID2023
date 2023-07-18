library(DSAIRM)
source("simulate_pkpdbacteria_ode.R")

# running new bacteria model with drug PK/PD to explore impact of drug
# simulation 1: drug present but not impact (Emax=0)
sim1 <- simulate_pkpdbacteria_ode(B = 10, I = 1, g = 1, Bmax = 1e+03, 
                            dB = 0.1, kI = 1e-3, r = 2e-3, dI = 1,
                            C0 = 1000, dC = 0.1, C50 = 500, k = 1, Emax = 0,
                            txstart = 10, txinterval = 10, tstart = 0, 
                            tfinal = 200, dt = 0.01)
sim1$yscale = "identity" #set to "log10" for log plot
generate_ggplot(list(sim1))


# simulation 1: drug present and working (Emax>0)
sim2 <- simulate_pkpdbacteria_ode(B = 10, I = 1, g = 1, Bmax = 1e+03, 
                                 dB = 0.1, kI = 1e-3, r = 2e-3, dI = 1,
                                 C0 = 1000, dC = 0.1, C50 = 500, k = 1, Emax = 1,
                                 txstart = 10, txinterval = 10, tstart = 0, 
                                 tfinal = 200, dt = 0.01)
sim2$yscale = "identity"
generate_ggplot(list(sim2))
