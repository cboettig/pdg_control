# file: training_prob2_collocation.R
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-27


alpha <- 1
K <- 1
r <- 1
T<- 1

# Discrete time state equation.  BH-style
# \[ x_{t+1} = x_t + f(t, x_t, u_t) \] 
f <- function(t, x, u){
  r * x ^ alpha / (1 + x ^ alpha / K) - x * u
}

# Define the utility function R
R <- function(t, x, u){
  gamma <- -2
  u ^ (1+gamma) / (1 + gamma)
}

# and the boundary condition cost \Phi
Phi <- function(x) (x != XT)*Inf 





# Define the Cost function
# \[ C(x_0, u(t)) = \Phi(x_T) + \sum_{t=0}^{T-1} R(t,x_t, u_t) \]
C <- function(x, u){
  Phi(x[length(x)]) + 
  sum(sapply(1:length(t), function(i) R(t[i], x[i], u[i])))
}

# The discrete time grid
#t <- seq(0, T, length.out = 10)

#U <- matrix(NA, ncol=2, nrow=length(t))
#J <- matrix(NA, ncol=length, nrow=length(t))

#J[


t = T-1.  x =



for(t in (T-1):0){
  for(x_t in x){

   <- R(t, x, u[i], 

  }
}


# Define the optimal cost-to-go:
#\[ J(t, x_t) = \min_{u_{t:T-1}}  C(x_t, u_{t:T-1} )  \]
J <- function



