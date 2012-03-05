<!--begin.rcode setup, echo=FALSE 
render_gfm()  
opts_knit$set(upload = TRUE)   
opts_knit$set(upload.fun = function(file){
   library(RWordPress) 
   uploadFile(file)$url
  })
## The real source code is externalized from this file:
end.rcode-->

<!--roptions dev="png", fig.width=7, fig.height=5, tidy=FALSE, warning=FALSE, message=FALSE, comment=NA, external=TRUE, cache=FALSE, cache.path="perfectpolicy/"-->

# perfect policy, imperfect implementation 
Compare to non-optimal, rule-of-thumb policy.

### Model setup 
Clear the workspace and load package dependencies: 
<!--begin.rcode libraries_, echo=FALSE
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
end.rcode-->

Define parameters
<!--begin.rcode parameters_
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition
end.rcode-->

Use log-normal noise functions
<!--begin.rcode noise_dists_
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
end.rcode-->

Chose the state equation / population dynamics function
<!--begin.rcode Myer
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point, for reference 
x0 <- K - sigma_g ^ 2 / 2 
end.rcode-->

<!--begin.rcode profit_fn
profit <- profit_harvest(price_fish = 1, 
                         cost_stock_effect = 0,
                         operating_cost = 0.1)
end.rcode-->

<!--begin.rcode grid
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
end.rcode-->


### The perfect policy 
Calculate the optimal policy
<!--begin.rcode precuationary 
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
end.rcode-->

### The imperfect implementation

Implementation the optimal polict with implementation noise 
<!--begin.rcode simulate_edited
sigma_i <- 0.4 
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
end.rcode-->

#### Outcome 
<!--begin.rcode tidy_
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
end.rcode-->

This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown (solid line), along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 
<!--begin.rcode fishstock_policy, fig.width=9
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
end.rcode-->

Calculate which crashed
<!--begin.rcode hascrashed
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
end.rcode-->
A total of <!--rinline sum(crashed$V1) --> crash.



### A non-optimal policy 
Let's adjust the optimal policy by a rule-of-thumb buffer, resulting in a non-optimal policy.
<!--begin.rcode safe_policy
buffer <- 0.1
safe_policy <- matrix(sapply(opt$D - buffer * length(h_grid), function(x) max(0, x)), ncol=dim(opt$D)[2])
end.rcode-->

This adds a <!--rinline 100*buffer--> % buffer below the optimal harvest rate. 


<!--begin.rcode simulate_edited_noisy
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,K,1), x_grid, h_grid, x0, safe_policy, z_g, z_m, z_i)
})
end.rcode-->

<!--begin.rcode ref.label="tidy_"
end.rcode-->

<!--begin.rcode fishstock_policy2, fig.width=9
policy <- melt(safe_policy)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
end.rcode-->

<!--begin.rcode ref.label="hascrashed"
end.rcode-->
A total of <!--rinline sum(crashed$V1) --> crash.



