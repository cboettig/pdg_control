## @knitr libraries
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)

## @knitr pars 
delta <- 0.05     # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 0       # bonus for satisfying the boundary condition




## @knitr noise_dists
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1




## @knitr BevHolt_
f <- BevHolt                # Select the state equation
pars <- c(1.5, 0.05)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference 
xT <- 0                     # boundary conditions
x0 <- K


## @knitr profit_
profit <- profit_harvest(price = 10, c0 = 30, c1 = 10)


## @knitr create_grid_
x_grid <- seq(0.01, 1.2 * K, length = gridsize)  
h_grid <- seq(0.01, 0.8 * K, length = gridsize)  


## @knitr fees
L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
asymmetric <- function(c2) function(h, h_prev)  c2 * max(h - h_prev, 0)
fixed <-  function(c2) function(h, h_prev) c2
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2


## @knitr loop over possibilities
penaltyfns <- list(L2=L2, L1=L1, asy=asymmetric, fixed=fixed)
policies <- 
lapply(penaltyfns, function(penalty){
  c2 <- seq(0.01, 3, length.out = 15)
  policies <- 
  lapply(c2, function(c2){
    SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                         profit, delta, reward=reward)
    policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                        profit, delta, reward, penalty = penalty(c2))
    sims <- lapply(1:100, function(i)
      simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D, profit=profit, penalty=penalty(c2))
      )
    dat <- melt(sims, id=names(sims[[1]])) 
    names(dat) <- c(names(sims[[1]]), "reps")
    dat
  })

  dat <- melt(policies, id = names(policies[[1]]))
  dat$L1 <- c2[dat$L1]
  names(dat) <- c(names(policies[[1]]), "c2")
 dat
})
dat <- melt(policies, id = names(policies[[1]]))
names(dat) <- c(names(policies[[1]]), "penalty")
c2 <-  seq(0.01, 3, length.out = 15)
c2_index <- sapply(dat$c2, function(c) which(c2 %in% c))
dat2 <- data.frame(dat, c2_index = c2_index)
dt <- data.table(dat2)
setkey(dt,penalty,c2_index,reps)


## Beautiful
vartrend <- 
sapply(c("L2", "L1", "fixed", "asy"), function(plty){
  sapply(1:15, function(i){
    pl <- quote(plty)
    mean(dt[penalty==eval(pl) & c2_index ==i,
         sd(harvest)/sd(fishstock), by=reps]$V1)
  })
})
df <- melt(data.frame(c2=c2, vartrend),id="c2") 
p1 <- ggplot(df) + geom_point(aes(c2, value, color=variable))+ geom_line(aes(c2, value, color=variable))


cortrend <- sapply(c("L2", "L1", "fixed", "asy"), function(plty){
  sapply(1:15, function(i){
    pl <- quote(plty)
    mean(dt[penalty==eval(pl) & c2_index ==i,
         cor(harvest, fishstock), by=reps]$V1)
  })
})
df <- melt(data.frame(c2=c2, cortrend),id="c2") 
p2 <- ggplot(df) + geom_point(aes(c2, value, color=variable)) + geom_line(aes(c2, value, color=variable))

discount <- (1-delta)^seq(1:50)

npv <- sapply(c("L2", "L1", "fixed", "asy"), function(plty){
  sapply(1:15, function(i){
    pl <- quote(plty)
    mean(dt[penalty==eval(pl) & c2_index ==i,
         sum(profit_fishing*discount), by=reps]$V1)
  })
})
df <- melt(data.frame(c2=c2, npv),id="c2") 
NPV0 <- mean(npv[1,])
p2 <- ggplot(df) + geom_point(aes(c2, (NPV0-value)/NPV0, color=variable)) + geom_line(aes(c2, (NPV0-value)/NPV0, color=variable))









save(list=ls(), file="quadratic.rda")
