###########################
# example code for Task 1 of the bacteria exploration app
# shows how to do it 2 different ways
###########################

# load package
library(DSAIRM)

###########################
# first approach using the model exploration app
# this app repeatedly calls the basic bacteria ode model
# this approach is easier to code, but not as flexible
###########################
# run app to loop/scan over parameter of interest
sim1 <- simulate_basicbacteria_modelexploration(B = 100, I = 10, 
                                                      g=2, Bmax=1e5, dB=1, 
                                                      k=1e-4, r=1e-4, dI=2, 
                                                      tstart = 0, tfinal = 20, dt = 0.1, 
                                                      samples = 10, parmin=2, parmax=10, 
                                                      samplepar='g',  pardist = 'lin')

# look at results returned from simulation
# only provides hard-coded outputs for peak and steady states, not full time series
str(sim1)
# pull out number of times steady state was not reached
N_not_steady = sum(!sim1$dat$steady)
print(N_not_steady)

###########################
# second approach
# instead of using the modelexploration helper function above
# this calls the basic bacteria simulation function directly
# it has a manual loop to run over all values of the parameter of interest
###########################

#store all results in a big list
# this is not really needed in general, but we can do it here so we have access to 
# the full time series information and can compute any outcome we want
all_results = list() 

# specify vector of parameter values to scan/loop over
g_vec = seq(2,10,length=10)

# loop over parameter of interest
for (n in 1:length(g_vec))
{
  # for each parameter value, run simulation model, get full time series back
  all_results[[n]] <- simulate_basicbacteria_ode(B = 100, I = 10, g = g_vec[n], 
                                                 Bmax = 1e+05, dB = 1, 
                                                 k = 1e-4, r = 1e-4, dI = 2, 
                                                 tstart = 0, tfinal = 20, dt = 0.1)
}
# look at large list containing all simulation results
str(all_results)

# post-process to get the quantities we are interested in
# this could be done inside the loop above
# I'm splitting it here to make it easier to follow each step

# to figure out if we reached steady state, we need to record bacteria at end
# and a few time steps before the end
Bfinal = rep(0,length(g_vec))
Bnotfinal = rep(0,length(g_vec))

# loop over simulations, 
# for each simulation pull out value of B at end and a few time steps before end
# here, we pull out the last 5 values of B for each simulation
# we record the final value and the one 5 time steps before the final
# if we reached steady state, they should be the same
for (n in 1:length(g_vec))
{
  Bfinal[n] = round(tail(all_results[[n]]$ts[,"B"],5)[5],0)
  Bnotfinal[n] = round(tail(all_results[[n]]$ts[,"B"],5)[1],0)
}

# this shows the values at the last time step and a few time steps prior
cat('Final values of B:', Bfinal)
cat('Not Final values of B:', Bnotfinal)

#number of times we didn't reach steady state 
#this is not exactly the same as above due to rounding errors
#when the time-series is returned, it is already rounded
sum(Bfinal-Bnotfinal != 0)

