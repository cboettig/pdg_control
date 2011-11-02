# file: training_prob2_collocation.U
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-27



# \[ x_{t+1} = x_t + f(x_t) - h_t \]
# F(x) = f(x) + x
#
# Discount rate \beta < 1
# \[ \max_{h_t} \sum_{t=0}^{\infty} \beta^t U(h_t)  \]
# s.t. \( x_{t+1} = F(x_t) -h_t \)
# x_0 = x0


# First-order conditions are:
#\[ U'(h_t) = \beta U'(h_{t+1}) F'(x_{t+1} ) \forall t \]

# \[ V(x) = \max_h U(h) + \beta V(x_{t+1} ) \]
# \[ V(x) = \max_h U(h) + \beta V( F(x) - h) \]



rm(list=ls()) # start clean

alpha <- 1
K <- 1
C <- 0.2
beta <- 0.9 # discounting rate (discrete)
gamma <- -2
T <- 4

# Discrete time state equation.  Ricker-style
# \[ x_{t+1} = f(t, x_t, h_t) \] 
f <- function(t, x, h){
  x - h
#  r * x ^ alpha / (1 + x ^ alpha / K) - x * h ## B-H style
#  x * exp(alpha * (1 - x / K) * (x - C) / K ) - h * x
}



p <- c(20,22,30,35)
c <- 0.05
# Define the utility function U
U <- function(t, x, h){
   p[t]*h - c*h^2 
#  h ^ (1+gamma) / (1 + gamma)
}

# and the boundary condition cost \Phi

#phi <- function(x) p*x - c*x^2 # harvest what remains 
phi <- function(x) 300

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
    h <- seq(0,1000,length=100)
    cost <- sapply(h, func)
    i <- which.max(cost)
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

y <- numeric(T)
y[1] <- 1000 # x_0 constraint/initial condition
h <- numeric(T)
for(t in 1:(T-1)){
  h[t] <- h_star(t, y[t])  
 y[t+1] = f(t,y[t], h[t])
}

png("optimal.png")
plot(1:T, y, pch=19, cex=1.5, ylim=c(0,1000) )
points(1:T, h, pch=18, col="red", cex=2.5)
legend("bottomleft", c("fish", "harvest"), pch=c(19,18), col=c("black", "red"))
dev.off()
#

