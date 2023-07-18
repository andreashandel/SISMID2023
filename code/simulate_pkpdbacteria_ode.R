simulate_pkpdbacteria_ode <- function(B = 100, I = 1, g = 1, Bmax = 1e+05, 
                                      dB = 0.5, kI = 1e-4, r = 1e-4, dI = 2,
                                   C0 = 1, dC = 1, C50 = 1, k = 1, Emax = 0,
                                   txstart = 10, txinterval = 1, tstart = 0, tfinal = 20, dt = 0.01)
{

  #function that specificies the ode model
  pkpdode <- function(t, y, parms)
  {
    with(
      as.list(c(y,parms)), #lets us access variables and parameters stored in y and parms by name
      {
        e=Emax*C^k/(C^k + C50) #drug efficacy, based on level of C
        dCdt = - dC*C
        dBdt = g*B*(1-B/Bmax) - dB*B - kI*B*I - e*B
        dIdt = r*B*I - dI*I
        list(c(dCdt, dBdt, dIdt))
      }
    ) #close with statement
  } #end function specifying the ODEs

  #function that specifies addition of drug at the indicated time
  adddrug <- function(t,y,parms)
  {
    y['C'] = y['C'] + parms['C0']
    return(y)
  }

  Y0 = c(C = 0, B = B, I = I);  #combine initial conditions into a vector - drug starts at zero in ODE
  timevec = seq(tstart, tfinal, by = dt); #vector of times for which solution is returned (not that internal timestep of the integrator is different)

  #combining parameters into a parameter vector
  pars = c(g = g, Bmax = Bmax, dB =dB, kI = kI, r=r, dI=dI, C0 = C0, dC=dC,C50 = C50, k = k, Emax = Emax);

  drugtimes = seq(txstart, tfinal, by = txinterval) #times at which drug is administered
  #this line runs the simulation, i.e. integrates the differential equations describing the infection process
  #the result is saved in the odeoutput matrix, with the 1st column the time, all other column the model variables
  #in the order they are passed into Y0 (which needs to agree with the order in virusode)
  odeoutput = deSolve::ode(y = Y0, times = timevec, func = pkpdode, events = list(func = adddrug, time = drugtimes), parms=pars, atol=1e-12, rtol=1e-12);

  #return result as list, with element ts containing the time-series
  result = list()
  result$ts = as.data.frame(odeoutput)
  return(result)
}
