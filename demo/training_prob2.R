# training_prob2.R

# population growth function: Discrete time Ricker with Allee effect
# Returns the function f as a closure
set_f <- function(gamma, C, K){ 
  function(x) x * exp( gamma*(1-x/K)*((x-C)/K))
}
# Initialize 
f <- set_f(1, 20, 100)

# Recursively apply f for t years
F <- function(N,t){ 
  i = 1
  n <- N
  iterate_t <- function(t){
    while(i<=t){ 
      n <- f(n)
      i<- i+1
    }
    n
  }
  sapply(t, iterate_t)
}

t <-1:20
plot(t, F(22,t))



Pi <- function(E_t, E_s, t){
# Define the initial probability distribution of population density
# SD comes at a given cost, sd = 1/E_s
  P_0 <- function(x) dnorm(x, mean=60, sd=1/E_s) 

  iterate_t <- function(t){ # in case t is a vector
    int <- function(N) E_t*N*F(N,t)*P_0(N)
    integrate(int, 0, 200)$value
  }
  sapply(t, iterate_t)
}


r <- 1 # discounting (interest) rate

# target function which we optimize
target <- function(pars){
  int <- function(t) Pi(pars["E_t"], pars["E_s"], t)*exp(-r*t)
  -integrate(int, 0, T)$value
}


E_s <- seq(.1, 2.1, by=.5)
E_t <- seq(.1, 2.1, by=.5)

out<-sapply(E_s, function(y){
  sapply(E_t, function(x){
    pars <- c(E_t=x,E_s=y)
    target(pars)
  })
})

# Whoops, solution is to harvest at maximum effort and spend nothing on sampling

#optim(pars, target)



