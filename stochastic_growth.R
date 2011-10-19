###############################################
#
# stochastic_growth.R: An R program to solve a simple stochastic
# growth model via eitehr value function or policy function 
# iteration, which is also known as Howard's Policy Improvement 
# Algorithm.
#
# Characteristics of the problem are:
#   (1) Infinite horizon
#   (2) Continuous state space and discrete time
#   (3) Income shock takes on two discrete states
#
# Assumptions are:
#   (1) Utility function: U(c) = c^(1-g)/(1-g)
#   (2) discount rate = 10%
#   (3) Consumption function: c = F(K) = (1/a)*K^a.  AKA the law of motion
#
# Based on code written by George Hall, July 2001
# 

function <- stochastic_growth(){
  value_iteration <-0
  if(value_iteration==1)
    print("Value Function Iteration")
  else 
    print("Policy Function Iteration")

# Economic parameter values

a <- 0.20     # production parameter for f(K)
b <- 0.10     # discount rate
b <- 1/(1+b)  # discount factor
g <- 1.5      # g>=1 in U(c) function



# Transition probability matrix
P <- matrix(c(.8,.2,.2,.8), nrow=2, byrow=T) # 
A_high <- 1.250 # high value for technology
A_low <- 0.50   # low value for technology

# Discretizing the continuous state variable, Capital
N <- 1501 # N is the grid-size, slower fro larger N

#########
# Typically, you would solve for the steady-state level and 
# use that to determine the range.  I know for this problem
# the steady-state is around K=1 in setting this range:
########
MaxK <- 1.5
MinK <- 0.1
K <- seq(MinK, MaxK, length.out = N) # uniform grid


#=========================================================#
# Tabulating consumption and utility for all combinations #
# K, K(t+1) for the grid                                  #
#=========================================================#
# R needs a meshgrid function. 
meshgrid <- function(a,b) {
    list(
           x=outer(b*0,a,FUN="+"),
           y=outer(b,a*0,FUN="+")
        )
} 
grid <- meshgrid(K, K)
Kt <- grid$x
Kt1 <- grid$y

# Consumption function all along the grid
consumeh <- A_high * (1/a) * Kt^a - Kt1
consumel <- A_low * (1/a) * Kt^a - Kt1

# Replace all negative consumption with NaN
# Not always necessary, e.g.g when u(c)=ln(c)
consumeh[consumeh < 0] <- NaN
consumel[consumel < 0] <- NaN

# Utility function at all grid points 
utilityh = consumeh ^ (1-g) / (1-g)
utilityl = consumel ^ (1-g) / (1-g)

# Replace the NaNs with -Inf to avoid choosing these values
utilityh[is.nan(utilityh)] <- -Inf
utilityl[is.nan(utilityl)] <- -Inf

# R doesn't have a repmat function; but easy to create:
repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

# Initializing iteration variables
iter <- 0
v <- matrix(0, nrow=N, ncol=2)
decis <- matrix(0, nrow=N, ncol=2)
VError<-10
Tol <- 1e-7


if(valueiteration==1){
#===============================================#
# Iteration on value function                   #
#===============================================#
  while(VError > Tol){
# Calculate at the value at high
    U <- utilityh+b*repmat(v %*% P[1,], 1, N)
    tvh <- apply(U, 2, max)
    tdecish <- apply(U, 2, which.max)
# Calculate at the value at low
    U <- utilityl+b*repmat(v %*% P[2,], 1, N)
    tvl <- apply(U, 2, max)
    tdecisl <- apply(U, 2, which.max)
# Combine
    tv <- c(tvh, tvl)
    tdecis <- c(tdecish, tdecisl)
# Check the value function
    Verror<- max(abs((tv-v)/v))
    v <- tv # Update the value function
    decis <- tdecis # track the index of decision
    iter <- iter + 1
    if(iter > 500) 
      break
  }
} else {
#===============================================#
# Solve for fixed point using policy iteration  #
#===============================================#
  while(VError > Tol){
# Calculate at the value at high
    U <- utilityh+b*repmat(v %*% P[1,], 1, N)
    tdecish <- apply(U, 2, which.max)
# Calculate at the value at low
    U <- utilityl+b*repmat(v %*% P[2,], 1, N)
    tdecisl <- apply(U, 2, which.max)
    tdecis <- c(tdecish, tdecisl)

# Solving out for the consumption that corresponds to the value 
# function using the index from above calculation;
# x <- argmax(U(x) + bPV)
 }
}



}


