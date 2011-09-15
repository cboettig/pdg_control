# training_prob2.R


N_true <- 50

# population growth function: Discrete time Ricker with Allee effect
# Returns the function f as a closure
set_f <- function(gamma, C, K, beta){ 
  function(x, E_h) x * exp( gamma*(1-x/K)*((x-C)/K)) - x*E_h*beta
}

# Initialize f.  Allee threshold 20, carrying capacity 100, .1 efficacy 
# (e.g. effort = 1/beta means catch all the stock)
f <- set_f(1, 5, 100, .1)



# Recursively apply f for t years
F <- function(N,t, E_h){ 
  i = 1
  n <- N
  iterate_t <- function(t){
    while(i<=t){ 
      n <- f(n, E_h)
      i<- i+1
    }
    n
  }
  sapply(t, iterate_t)
}


info_cost <- 1


Pi <- function(E_h, E_s, t){
# Define the initial probability distribution of population density
# Increased sampling effort reduces the the variance 
  P_0 <- function(x) dlnorm(x, mean=log(N_true), sd=1/E_s)

  iterate_t <- function(t){ # in case t is a vector
    int <- function(N) E_h*N*F(N,t,E_h)*P_0(N)
    sum(sapply(0:200, int))
#    integrate(int, 0, 200)$value
  }

  info_cost_per_time <- function(t){
# pay only on first time interval
    out<- numeric(length(t))
    out[1] <- E_s*info_cost
    out
  }

  sapply(t, iterate_t) - info_cost_per_time(t)
}


r <- .5 # discounting (interest) rate
T <- 40 # time horizon

# target function which we optimize
target <- function(pars){
  int <- function(t) Pi(pars["E_h"], pars["E_s"], t)*exp(-r*t)
  -sum(sapply(1:T, int))
#  -integrate(int, 0, T)$value
}

E_s <- seq(1, 20,length=6)
E_h <- seq(1, 10, length=10)
## E_h rows, E_s cols
gridsearch <- function(E_s, E_h){
  out<-sapply(E_s, function(y){
    sapply(E_h, function(x){
      pars <- c(E_h=x,E_s=y)
      target(pars)
    })
  })
  out
}



pars <- c(E_h=.5, E_s=.5)
o <- optim(pars, target)

#optimal solution is:
o$par

par(mfrow=c(1,2)) # space for two side-by-side plots
## what does the optimal strategy do to the fish population?
plot(1:T, F(N_true,t,o$par["E_h"]), main="Optimial fish stock over time")
## What does the optimal starting information look like:
P_0 <- function(x) dlnorm(x, mean=log(N_true), sd=1/o$par["E_s"])
curve(P_0, 0, 120, main="Optimal Sampling")
## enough to place most of the weight above the allee threshold 



