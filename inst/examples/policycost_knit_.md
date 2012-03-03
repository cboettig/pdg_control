<!--roptions dev='png', fig.width=7, fig.height=5, tidy=FALSE, warning=FALSE, message=FALSE, comment=NA, cache.path="policycost/", cache=FALSE-->

<!--begin.rcode setup, include=FALSE
render_gfm()  
opts_knit$set(upload = TRUE)   
opts_knit$set(upload.fun = function(file){
   library(RWordPress) 
   uploadFile(file)$url
  })
read_chunk("Reed.R")
end.rcode-->


# Policy Costs 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

<!--begin.rcode libraries, echo=FALSE
end.rcode-->


This example illustrates the impact of adding a cost to changing the harvest level between years 

### Define all parameters 
<!--begin.rcode parameters
delta <- 0.01     # SMALLER economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition
end.rcode-->

we'll use log normal noise functions
<!--begin.rcode noise_dists
end.rcode-->


Chose the state equation / population dynamics function
<!--begin.rcode Myer
end.rcode-->

Our initial condition is the equilibrium size (note the stochastic deflation of mean)
<!--begin.rcode initx
end.rcode-->

and we use a harvest-based profit function with default parameters
<!--begin.rcode profit_
profit <- profit_harvest(price_fish = 2, cost_stock_effect = 0,
 operating_cost = 0.1)
end.rcode-->

Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 
<!--begin.rcode create_grid
end.rcode-->

### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 
<!--begin.rcode determine_SDP_matrix
end.rcode-->
### Find the optimum by dynamic programming 
We use Bellman's dynamic programming algorithm to compute the optimal solution for all possible trajectories, ignoring potential policy costs as before.  We will later use this solution to compare against the optimal solution with policy costs.
<!--begin.rcode find_dp_optim 
end.rcode-->

A modified algorithm lets us include a penalty of magnitude `P` and a functional form that can be an `L1` norm, `L2`  norm, `asymmetric` L1 norm (costly to lower harvest rates), fixed cost, or `none` (no cost).  Here is an asymmetric norm example.  Note that this calculation is considerably slower. 
<!--begin.rcode policycost_optim
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, P = .5, penalty = "asym")
end.rcode-->


### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`

<!--begin.rcode simulate_policy
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D)
  )
end.rcode-->


## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.
<!--begin.rcode tidy
end.rcode-->

### Plots 

Compare the optimal policy that involves this cost:
<!--begin.rcode policy_cost_vis, fig.width=10
policy <- melt(policycost$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=xT, slope=0, lty=2)
p5 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1, data=dt)
end.rcode-->

Against the policy with no cost: 
<!--begin.rcode no_policy_cost_vis, fig.width=10
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$alternate) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)  
p6 + geom_line(aes(time, alternate, group = reps), alpha = 0.1, data=dt)
end.rcode-->


Compare dynamics on a single replicate to see how this policy differs from the Reed policy. 
<!--begin.rcode plot_rep2
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 

end.rcode-->

## Alternate policy cost models 

#### L2 norm
<!--begin.rcode policycost_optim_l2
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, P = .5, penalty = "L2")
end.rcode-->
<!--begin.rcode policy_cost_vis_l2, fig.width=10
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D)
)

dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice


policy <- melt(policycost$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=xT, slope=0, lty=2) +
  geom_line(aes(time, alternate, group = reps), alpha = 0.1, data=dt)
end.rcode-->


#### L1 norm
<!--begin.rcode policycost_optim_l1
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, P = .5, penalty = "L1")
end.rcode-->
<!--begin.rcode policy_cost_vis_l1, fig.width=10
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D)
)

dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice



policy <- melt(policycost$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=xT, slope=0, lty=2) +
  geom_line(aes(time, alternate, group = reps), alpha = 0.1, data=dt)
end.rcode-->


