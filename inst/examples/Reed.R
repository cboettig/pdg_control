
## @knitr libraries
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)

## @knitr parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition

## @knitr noise_dists
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1



## @knitr BevHolt 
f <- BevHolt                # Select the state equation
pars <- c(2, 4)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity 
xT <- 0                     # boundary conditions

## @knitr Myer
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point, for reference 

## @knitr RickerAllee
f <- RickerAllee
K <- 4 
xT <- 1 # final value, also allee threshold
pars <- c(1, K, xT) 


## @knitr initx
x0 <- K - sigma_g ^ 2 / 2 

## @knitr profit 
profit <- profit_harvest(price_fish = 1, cost_stock_effect = 0,
 operating_cost = 0.1 * price)


## @knitr create_grid
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  

## @knitr determine_SDP_matrix
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )

## @knitr stochastic
require(snowfall) 
sfInit(parallel=TRUE, cpu=4)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)

## @knitr find_dp_optim 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)

## @knitr simulate
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})

## @knitr tidy
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice

## @knitr plot_rep
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 

## @knitr fishstock 
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)


## @knitr harvest
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")

## @knitr escapement
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")

## @knitr crashed
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]

## @knitr profits
dt <- data.table(dt, id=1:dim(dt)[1])
profits <- dt[, profit(fishstock, harvest), by=id]

## @knitr join
setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
setkey(dt, reps)

## @knitr total_profit
total_profit <- dt[,sum(profits), by=reps]
total_profit <- total_profit + rewarded$V1 * reward 

## @knitr joinmore
setkey(total_profit, reps)
setkey(crashed, reps)
setkey(rewarded, reps)
dt <- dt[total_profit]
dt <- dt[crashed]
dt <- dt[rewarded]
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))

## @knitr profit_by_time
stats <- dt[ , mean_sdl(profits), by = time]
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)

## @knitr totals
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)

## @knitr quantile
quantile_me <- function(x, ...){
  q <- quantile(x, ...)
  class <- character(length(x))
  for(i in 1:length(q))
    class[x > q[i] ] <- i
  class
}
q <- data.table(reps=total_profit$reps, quantile=quantile_me(total_profit$V1))
setkey(q, reps)
dt <- dt[q]

## @knitr winners_losers
ggplot(subset(dt, quantile %in% c(1,4))) + 
  geom_line(aes(time, fishstock, group = reps, color=quantile), alpha = 0.6) 

## @knitr policyvis
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p5


## @knitr policyvis2
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1, data=dt)

