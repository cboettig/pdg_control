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

We can also look at the harvest dynamics:
<!--begin.rcode harvest
end.rcode-->

This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement and see if it is less variable than fish stock, and if it is near Reed's S: 
<!--begin.rcode escapement
end.rcode-->

### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?
<!--begin.rcode crashed
end.rcode-->

A total of 
<!--begin.rcode crashedsum
sum(crashed$V1)
end.rcode--> 
crash.

Let's compute the profits at each time-step for each replicate. 
Using `data.table` to evaluate our profit function over the stock and harvest levels requires indexing our data:

<!--begin.rcode profits
end.rcode-->

Merging this calculation back into our data table using fast join (needs to define 'id' as a key on which to match things up though). 
<!--begin.rcode join
end.rcode-->

Compute total profit by summing over each timeseries (including the reward for satisfying the terminal boundary condition, if any). 

<!--begin.rcode total_profit
end.rcode-->


Add these three columns to the data.table (fast join and re-label):
<!--begin.rcode joinmore
end.rcode-->



#### Profit plots
Since the optimal strategy maximizes expected profit, it may be more useful to look at the distribution statistics of profit over time:
<!--begin.rcode profit_by_time
end.rcode-->


Total profits
<!--begin.rcode totals
end.rcode-->

## Compare to a non-optimal solution
Compare another model, that likewise assumes no implementation error, and also makes a mistake in its estimate of the growth parameter, making it conservative rather than optimal.

<!--begin.rcode model2_guesses
pars <- c(1, K, 1.5)
end.rcode-->


<!--begin.rcode find_dp_optim 
end.rcode-->

For the simulated implementation, we add the same implementation error back, and we restore biological growth noise to it's true value
<!--begin.rcode implement_error_again
sigma_i <- 0.4
end.rcode-->


### Simulate 
Now we simulate as before
<!--begin.rcode simagain, ref.label="simulate"
end.rcode-->

## Summarize and plot the results                                                  
Using the code above, recreate the plots for this policy and simulation: 
<!--begin.rcode tidy2, ref.label="tidy"
end.rcode-->

### Plots 
<!--begin.rcode ref.label="fishstock"
end.rcode-->

<!--begin.rcode ref.label="harvest"
end.rcode-->

<!--begin.rcode ref.label="escapement"
end.rcode-->

### Computing additional statistics about the data
<!--begin.rcode ref.label="crashed"
end.rcode-->
A total of 
<!--begin.rcode crashedsum
sum(crashed$V1)
end.rcode--> 
crash.

<!--begin.rcode ref.label="profits"
end.rcode-->

<!--begin.rcode ref.label="join"
end.rcode-->

<!--begin.rcode ref.label="total_profit"
end.rcode-->

<!--begin.rcode ref.label="joinmore"
end.rcode-->

#### Profit plots
<!--begin.rcode ref.label="profit_by_time" 
end.rcode-->

<!--begin.rcode ref.label="totals" 
end.rcode-->


