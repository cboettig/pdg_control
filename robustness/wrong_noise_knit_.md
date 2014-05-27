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

<!--roptions dev="png", fig.width=7, fig.height=5, tidy=FALSE, warning=FALSE, comment=NA, external=TRUE, cache=TRUE, cache.path="wrong_noise/"-->

# Reed Model, when in reality growth noise is slightly larger
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in:
 
   Sethi, G., Costello, C., Fisher, A., Hanemann, M., and Karp, L. (2005). 
   Fishery management under multiple uncertainty. Journal of Environmental
   Economics and Management, 50(2), 300-318. doi:10.1016/j.jeem.2004.11.005

   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
   and Deterministic Harvesting Models. Journal of Environmental 
   Economics and Management. 6: 350-363.

 
  Fish population dynamics:
 \\( X_{t+1} = Z_n f(X_n) \\)


Clear the workspace and load package dependencies: 
<!--begin.rcode libraries, echo=FALSE
end.rcode-->




### Define all parameters 
<!--begin.rcode parameters
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
<!--begin.rcode profit
end.rcode-->


Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 

<!--begin.rcode create_grid
end.rcode-->


### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. 

<!--begin.rcode determine_SDP_matrix
end.rcode-->


<!--begin.rcode stochastic, eval=FALSE, include=FALSE
end.rcode-->

### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.
<!--begin.rcode find_dp_optim 
end.rcode-->

### Reality is just a bit more noisy than we think

<!--begin.rcode reality
sigma_g <- 1.2 * sigma_g
end.rcode-->


### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.
<!--begin.rcode simulate
end.rcode-->




## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.
<!--begin.rcode tidy
end.rcode-->

### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above. 
<!--begin.rcode plot_rep
end.rcode-->


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


#### Add discrete classes by total profit

Sometimes I'd like to color code the replicates by profit, to see if there are particular patterns in stock dynamics of the most profitable and least profitable lines.  Adding discrete profit classes to the data table makes this possible:
<!--begin.rcode quantile
end.rcode-->

Then we can plot the fishstock trajectories, indicating which derive the highest and smallest profits by color code: 
<!--begin.rcode winners_losers
end.rcode-->

### Visualizing the optimal policy
<!--begin.rcode policyvis, fig.width=9
end.rcode-->

<!--begin.rcode policyvis2, fig.width=10
end.rcode-->

