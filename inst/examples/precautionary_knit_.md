<!--begin.rcode setup, echo=FALSE 
render_gfm()  
opts_knit$set(upload = TRUE)   
opts_knit$set(upload.fun = function(file){
   library(RWordPress) 
   uploadFile(file)$url
  })
## The real source code is externalized from this file:
read_chunk("Reed.R")
end.rcode-->

<!--roptions dev="png", fig.width=7, fig.height=5, tidy=FALSE, warning=FALSE, message=FALSE, comment=NA, external=TRUE, cache=FALSE, cache.path="perfectpolicy/"-->

# Precautionary Principle


## Model setup 

Clear the workspace and load package dependencies: 
<!--begin.rcode libraries, echo=FALSE
end.rcode-->

Define parameters
<!--begin.rcode parameters
end.rcode-->

Use log-normal noise functions
<!--begin.rcode noise_dists
end.rcode-->

Chose the state equation / population dynamics function
<!--begin.rcode RickerAllee
end.rcode-->

Our initial condition is the equilibrium size (note the stochastic deflation of mean)
<!--begin.rcode initx
end.rcode-->

and we use a harvest-based profit function with default parameters
<!--begin.rcode profit
end.rcode-->

Set up discrete grids for stock size and havest levels (which will use same resolution as for stock). 
<!--begin.rcode create_grid
end.rcode-->


### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. 

Here we assume a precautionary high value of the threshold parameter.
<!--begin.rcode determine_SDP_matrix_edited
  SDP_Mat <- determine_SDP_matrix(f, c(1,4,2), x_grid, h_grid, sigma_g )
end.rcode-->

### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.
<!--begin.rcode find_dp_optim 
end.rcode-->

### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above, with a threshold value that is actually lower than we assumed. 
<!--begin.rcode simulate_edited
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,K,1), x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
end.rcode-->

## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.
<!--begin.rcode tidy
end.rcode-->

### Plots 
This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 
<!--begin.rcode fishstock 
end.rcode-->

<!--begin.rcode crashed
end.rcode-->
A total of <!--rinline sum(crashed$V1) --> crash.



### So is this more robust?
Here's the policy under higher intrinsic noise than it was designed for:
<!--begin.rcode simulate_edited_noisy
sigma_g <- .25
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,K,1), x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
end.rcode-->

Make data tidy (melt), fast (data.tables), and nicely labeled.
<!--begin.rcode ref.label="tidy"
end.rcode-->

<!--begin.rcode ref.label="fishstock"
end.rcode-->

<!--begin.rcode ref.label="crashed"
end.rcode-->
A total of <!--rinline sum(crashed$V1) --> crash.


## Compare to the optimal policy performance
<!--begin.rcode
pars <- c(1,K,1)
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g=.2)
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
sigma_g <- .25
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
crashed <- dt[time==as.integer(OptTime-1), fishstock == 0, by=reps]
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
end.rcode-->

A total of <!--rinline sum(crashed$V1) --> crash.

