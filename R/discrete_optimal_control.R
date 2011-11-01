# file: training_prob2_collocation.U
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-27



# \[ k_{t+1} = k_t + f(k_t) - h_t \]
# F(x) = f(k) + k
#
# Discount rate \beta < 1
# \[ \max_{h_t} \sum_{t=0}^{\infty} \beta^t U(h_t)  \]
# s.t. \( k_{t+1} = F(k_t) -h_t \)
# k_0 = k0


# First-order conditions are:
#\[ U'(h_t) = \beta U'(h_{t+1}) F'(k_{t+1} ) \forall t \]

# \[ V(k) = \max_h U(h) + \beta V(k_{t+1} ) \]
# \[ V(k) = \max_h U(h) + \beta V( F(k) - h) \]



rm(list=ls()) # start clean

alpha <- 1
K <- 1
C <- 0.2
beta <- 0.95 # discounting rate (discrete)
gamma <- -2
T <- 10

# Discrete time state equation.  Ricker-style
# \[ x_{t+1} = f(t, x_t, h_t) \] 
f <- function(t, x, h){
  if(is.na(h)) 
    recover()
#  r * x ^ alpha / (1 + x ^ alpha / K) - x * h ## B-H style
  x * exp(alpha * (1 - x / K) * (x - C) / K ) - h * x
}




# Define the utility function U
U <- function(t, x, h){
  h ^ (1+gamma) / (1 + gamma)
}

# and the boundary condition cost \Phi
phi <- function(x) if(x > 0.2) 0 else 1e6

Phi <- function(x) {
  sapply(x, 
    function(X){
      if(is.na(X))
        out <- Inf
      else if(X > .2) 
        out <- 0
      else
        out <- Inf
      out
    })
}

J <- function(t,x){
  if(t < T)
    out <- U(t,x, h_star(t,x)) + beta*J(t+1, f(t,x, h_star(t,x)))
  else 
    out <- phi(x)
  out
}
  
h_star <- function(t,x){
    func <- function(h) U(t, x, h) + beta*J(t+1, f(t,x,h))
#    optimize(f=func, interval=c(0,1))[[1]]
    h <- seq(0,1,length=10)
    cost <- sapply(h, func)
    i <- which.max(cost)
    if(is.na(h[i])){
      print(paste("failed i=", i, "length(cost) =", length(cost)))
      }
    h[i]
}

#
#### Constant fishing effort check  ###
#y <- numeric(T)
#y[1] <- K # constraint/initial condition
#for(t in 1:(T-1))
# y[t+1] = f(t,y[t], .2)
#plot(1:T, y, pch=19)
#
#
### Optimal Control Solution ##
#
#y <- numeric(T)
#y[1] <- K # x_0 constraint/initial condition
#h <- numeric(T)
#for(t in 1:(T-1)){
#  h[t] <- h_star(t, y[t])  
# y[t+1] = f(t,x[t], h[t])
#}
#plot(1:T, y, pch=19)
#points(1:T, h, pch=19, col="red")
#
#

