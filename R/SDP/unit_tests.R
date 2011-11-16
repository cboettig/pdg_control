# test the performance of SDP_by_simulation 

source("stochastic_dynamic_programming.R")

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
sigma_g <- 0.2      # Noise in population growth
gridsize <- 100   # gridsize (discretized population)
sigma_m <- .2  # 
sigma_i <- .0 # 
interval <- 1

# Chose the state equation / population dynamics function
f <- BevHolt 
pars <- c(2,4)  
K <- (pars[1]-1)/pars[2] 
xT <- 0 
e_star <- 0

x_grid <- seq(0, 2*K, length=gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same res as stock

SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)

require(snowfall)
sfInit(parallel=TRUE, cpu=16)
sim_SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, sigma_g, 0, 0)

x <- matrix(NA, length(x_grid), 10) 
x[,1] <- 0
x0 <- .2
x[which.min(abs(x_grid - x0)),1] <- 100
y <- x
for(i in 1:9){
  x[,i+1] <- t(sim_SDP_Mat[[1]]) %*% x[,i]
  y[,i+1] <- t(SDP_Mat[[1]]) %*% y[,i]
}

save(list=ls(), file = "unit_tests.rda")

png("test1.png")
  barplot(x[,2], col=rgb(0,0,1,.5))
  barplot(y[,2], add=T, col=rgb(1,0,0,.5))
dev.off()
png("test2.png")
  barplot(x[,5], col=rgb(0,0,1,.5))
  barplot(y[,5], add=T, col=rgb(1,0,0,.5))
dev.off()

require(socialR)
upload("test*.png", script="unit_tests.R", tag="PDG_Control")

