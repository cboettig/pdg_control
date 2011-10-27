# file: training_prob2_collocation.R
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-24

# Libraries
require(bvpSolve)


#######################################################################
# Set up model and parameters                                         #
#######################################################################

## Biological Parameters (state equation)
alpha <- 1  # Population growth rate
K <- 1      # carrying capacity
C <- 0.2    # allee theshold

## Economic Parameters 
gamma <- -2 #  
rho <- 0.5 # discounting

## Store these in a nice vector.  Won't name vars, as is slower
pars <- c(alpha, K, C, gamma, rho)

## Boundary Conditions
X0 <- K     # Starting population size
XT <- C     # Ending population size


# y is a vector of (x,h)', the fish population and harvest level

## State Equation(s): Fish population dynamics
f <- function(t, y, p){
  # Rename explicitly so equation is easier to read but still fast.
  x <- y[1]
  h <- y[2]
  alpha <- pars[1]
  K <- pars[2]
  C <- pars[3]
  x * alpha * ((K - x) / K) * ((x - C) / K) - h*x
}
## Derivative with respect to state x
df <- function(t, y, p){
  x <- y[1]
  h <- y[2]
  alpha <- pars[1]
  K <- pars[2]
  C <- pars[3]
  - alpha * ( C* K - 2 * x * K - 2 * x * C + 3 * x ^ 2) / K ^ 2 - h
}

######################################################################
# Solve the ODE system of fish dynamics at fixed harvest level       #
######################################################################
fish <- function(t,y,p){
  dy1 <- f(t, y, p) 
  dy2 <- 0                # constant harvest level 
  list(c(dy1, dy2), NULL)
}
jac <- function(t,y,p){
  matrix(c(df(t, y, p), 0,
           0,           0),
         2,2, byrow=T)
}
y0 <- c(0.5, .2)
t <- seq(0, 10, by=0.1)
out <- lsoda(y=y0, times=t, func=fish, parms=p, rtol=1e-4, jacfun=jac)
plot(out[,1], out[,2], type="l", xlab="time", ylab="fish stock", lwd=3)


#######################################################################
# Specify & solve the Boundary Value Problem                          #
#######################################################################

#' fun defines the bvp system to be solved
#' @param t time variable 
#' @param y a vector of the system state y[1],y[2] = (x,h)
#' @param p parameters: c(alpha, K, C, gamma, rho)
fun <- function(t,y,p){
  gamma <- p[4]
  rho <- p[5]

  dy1 <- f(t, y, p) 
  dy2 <- (1 / gamma) * (rho - df(t,y,p))

  list(c(dy1, dy2))
}

## Solve by shooting  ##
t <- seq(0, 10, by=.1)
sol1 <- bvpshoot(yini = c(X0,NA), yend = c(XT, NA),
                 x = t, func = fun, guess = 0, parms=pars)

## Solve by collocation ##
## (using the shooting solution as the guess)
sol3 <- bvpcol(yini = c(X0, NA), yend = c(XT, NA),
               x = t, func = fun, parms=pars, 
               xguess = sol1[,1], yguess=t(sol1[,-1]))


plot(sol3[,1], sol3[,2], type="l", lwd=3, ylim=c(0,1), xlab="time", ylab="fraction of K")
lines(sol3[,1], sol3[,3], col="red", lwd=3)
legend("topright", c("fish stock", "harvest effort"), lty=1, col=c("black", "red"))

## check solution quality 
diagnostics(sol3)



## save output to file
png("collocation.png")
  plot(sol3[,1], sol3[,2], type="l", lwd=3, ylim=c(0,1), xlab="time", ylab="fraction of K")
  lines(sol3[,1], sol3[,3], col="red", lwd=3)
  legend("topright", c("fish stock", "harvest effort"), lty=1, col=c("black", "red"))
dev.off()
# sol2 <- bvptwp(yini = c(1, NA), yend = c(.35, NA),
#               x = x, func = fun, parms=pars)



