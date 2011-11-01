# file: training_prob2_collocation.R
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-24

# Libraries
require(bvpSolve)

# load the functions
source("training_prob2_collocation.R")

##################################################################
# Set up model and parameters                                    #
##################################################################


## Biological Parameters (state equation)
# alpha - population growth rate
# K - carrying capacity
# C - allee theshold
## Economic Parameters 
# gamma - utility function (diminishing returns) 
# rho - economic discounting
pars <- c(alpha=1, K=1, C=0.2, gamma=-2, rho=0.5)


## Boundary Conditions
X0 <- K     # Starting population size
XT <- C     # Ending population size

y0 <- c(0.5*K, 0) # 
t <- seq(0, 10, by=0.1)
out <- lsoda(y=y0, times=t, func=fish, parms=p, rtol=1e-4, jacfun=jac)
plot(out[,1], out[,2], type="l", xlab="t", ylab="fish stock", lwd=3)


## Solve by shooting  ##
t <- seq(0, 10, by=.1)
sol1 <- bvpshoot(yini = c(X0,NA), yend = c(XT, NA),
                 x = t, func = fun, guess = 0, parms=pars)

## Solve by collocation ##
## (using the shooting solution as the guess)
sol3 <- bvpcol(yini = c(X0, NA), yend = c(XT, NA),
               x = t, func = fun, parms=pars, 
               xguess = sol1[,1], yguess=t(sol1[,-1]))

## check solution quality 
diagnostics(sol3)

## save output to file
png("collocation.png")
  time <- sol3[,1]
  fishstock <- sol3[,2] 
  harvesteff <- sol3[,3] 

  ylim <- c(min(fishstock, harvesteff), max(fishstock, harvesteff)) 
  plot(time, fishstock, type="l", lwd=3, 
       ylim=ylim, xlab="time", ylab="fraction of K")
  lines(sol3[,1], harvesteff, col="red", lwd=3)
  legend("topright", c("fish stock", "harvest effort"), lty=1, 
         col=c("black", "red"))

mtext(expression(bold("Optimal Harvest for:")), 
       line=2.5,cex=1.15)
mtext(bquote(paste(alpha==.(pars["alpha"]),", ", 
                   K==.(pars["K"]),", ", 
                   C==.(pars["C"]),", ", 
                   gamma==.(pars["gamma"]),", ", 
                   rho==.(pars["rho"]), sep="")), 
      line=1.25,cex=1.15) 
dev.off()
# sol2 <- bvptwp(yini = c(1, NA), yend = c(.35, NA),
#               x = x, func = fun, parms=pars)


