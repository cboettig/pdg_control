






# Asymmetric Policy Costs 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

## Setup the system



```r
source("setup.R")
```




### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 



```r
    SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




### Find the optimum by dynamic programming 

I've updated the algorithm to allow an arbitrary penalty function. Must be a function of the harvest and previous harvest. 



```r
c2 <- 4
penalty <- asymmetric(c2)

policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, penalty = penalty)
```





### Simulate 

Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`



```r
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D, profit=profit, penalty=penalty)
  )
```




Make data tidy (melt), fast (data.tables), and nicely labeled.



```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




### Plots 

A single replicate, alternate dynamics should show the Reed optimum, while harvest/fishstock should show the impact of having policy costs. 



```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, alternate)) +
  geom_line(aes(time, fishstock), col="darkblue") +
  geom_line(aes(time, harvest), col="purple") + 
  geom_line(aes(time, harvest_alt), col="darkgreen") 
```

![plot of chunk rep1](http://farm8.staticflickr.com/7115/7134078653_c040e0a5b8_o.png) 


A second replicate



```r
ggplot(subset(dt,reps==2)) +
  geom_line(aes(time, alternate)) +
  geom_line(aes(time, fishstock), col="darkblue") +
  geom_line(aes(time, harvest), col="purple") + 
  geom_line(aes(time, harvest_alt), col="darkgreen") 
```

![plot of chunk rep2](http://farm8.staticflickr.com/7066/6987994734_5c8b354511_o.png) 


We can visualize the equilibrium policy for each possible harvest:



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][,1])
ggplot(melt(policy)) + 
  geom_point(aes(h_grid[Var2], (x_grid[Var1]), col=h_grid[value]-h_grid[Var2])) + 
    labs(x = "prev harvest", y = "fishstock") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk policy](http://farm8.staticflickr.com/7176/7134079111_02cd2cfa40_o.png) 


Here we plot previous harvest against the recommended harvest, coloring by stocksize.  Note this swaps the y axis from above with the color density.  Hence each x-axis value has all possible colors, but they map down onto a subset of optimal harvest values (depending on their stock). 



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][,1])
ggplot(melt(policy)) + 
  geom_point(aes(h_grid[Var2], (h_grid[value]), col = x_grid[Var1]), position=position_jitter(w=.005,h=.005), alpha=.5) + 
    labs(x = "prev harvest", y = "harvest") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk harvestchanges](http://farm9.staticflickr.com/8026/6987995174_88a0edb1b4_o.png) 


## Profits



```r
save(list=ls(), file="stochastic_norms.rda")
```



