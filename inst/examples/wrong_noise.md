




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






### Define all parameters 



we'll use log normal noise functions




Chose the state equation / population dynamics function





Our initial condition is the equilibrium size (note the stochastic deflation of mean)



and we use a harvest-based profit function with default parameters




Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 





### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. 








### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.



### Reality is just a bit more noisy than we think



```r
sigma_g <- 1.2 * sigma_g
```



```
Error: object 'sigma_g' not found
```





### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.






## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.



### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above. 




This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 



We can also look at the harvest dynamics:



This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement and see if it is less variable than fish stock, and if it is near Reed's S: 





### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?



Let's compute the profits at each time-step for each replicate. 
Using `data.table` to evaluate our profit function over the stock and harvest levels requires indexing our data:




Merging this calculation back into our data table using fast join (needs to define 'id' as a key on which to match things up though). 



Compute total profit by summing over each timeseries (including the reward for satisfying the terminal boundary condition, if any). 





Add these three columns to the data.table (fast join and re-label):





#### Profit plots
Since the optimal strategy maximizes expected profit, it may be more useful to look at the distribution statistics of profit over time:




Total profits




#### Add discrete classes by total profit

Sometimes I'd like to color code the replicates by profit, to see if there are particular patterns in stock dynamics of the most profitable and least profitable lines.  Adding discrete profit classes to the data table makes this possible:



Then we can plot the fishstock trajectories, indicating which derive the highest and smallest profits by color code: 





