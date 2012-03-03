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

# Perfect Policy, Imperfect Implementation
In this example we will compare against an imperfect policy that has overestimated the noise, rather than the threshold position.  

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
<!--begin.rcode determine_SDP_matrix
end.rcode-->

### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.
<!--begin.rcode find_dp_optim 
end.rcode-->

### The optimal policy is implemented imperfectly
We add implementation noise: an imperfect implementation (though a symmetric one -- on average the implementation is not worse than assumed by the optimal solution, it is simply variable). 
<!--begin.rcode implementation_errors
sigma_i <- 0.4
end.rcode-->

### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above, but with this additional implementation error
<!--begin.rcode simulate
end.rcode-->

## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.
<!--begin.rcode tidy
end.rcode-->

### Plots 
This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 
<!--begin.rcode fishstock 
end.rcode-->

### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?
<!--begin.rcode crashed
end.rcode-->

A total of <!--rinline sum(crashed$V1) --> crash.



## Compare to a non-optimal solution
Compare another model, that likewise assumes no implementation error, and also makes a mistake in its estimate of the growth parameter, making it conservative rather than optimal.


<!--begin.rcode redoSDP
sigma_i <- 0
sigma_g <- .4 
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
end.rcode-->

<!--begin.rcode redoOpt
nonopt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
end.rcode-->


### Simulate 
For the simulated implementation, we add the same implementation error back, and we restore biological allee threshold to it's true value. 
<!--begin.rcode simagain
sigma_i <- .4
sigma_g <- .2
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, nonopt$D, z_g, z_m, z_i)
})
end.rcode-->

## Summarize and plot the results                                                  
Using the code above, recreate the plots for this policy and simulation: 
<!--begin.rcode tidy2, ref.label="tidy"
end.rcode-->

### Plots 
<!--begin.rcode ref.label="fishstock"
end.rcode-->

### Computing additional statistics about the data
<!--begin.rcode ref.label="crashed"
end.rcode-->
A total of <!--rinline sum(crashed$V1) --> crash.


