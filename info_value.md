---
layout: page
---





{% highlight text %}
Error: there is no package called 'knitcitations'
{% endhighlight %}




# Calculating the value of information
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in 

{% highlight text %}

Error in loadNamespace(name) : there is no package called 'knitcitations'

{% endhighlight %}

.
 Compute the optimal solution under different forms of uncertainty and compare the results.  





Define noise parameters 


Chose the state equation / population dynamics function



{% highlight r %}
f <- BevHolt
{% endhighlight %}




Note that the `pdg_control` pacakge already has a definition for the `BevHolt` function, (typing the function name prints the function)



{% highlight r %}
BevHolt
{% endhighlight %}



{% highlight text %}
function (x, h, p) 
{
    x <- max(0, x - h)
    A <- p[1]
    B <- p[2]
    sapply(x, function(x) {
        x <- max(0, x)
        max(0, A * x/(1 + B * x))
    })
}
<environment: namespace:pdgControl>
{% endhighlight %}




That is, \\( f(x,h) = \frac{A x}{1 + B x} \\)

Of course we could pass in any custom function of stocksize `x`, harvest `h` and parameter vector `p` in place of `BevHolt`.  Note that we would need to write this function explicitly so that it can take vector values of `x` (i.e. uses `sapply`), an annoying feature of `R` for users comming from Matlab.  


We must now define parameters for the function.  Note that the positive stationary root of the model is given by \\( \frac{A-1}{B} \\), which we'll store for future reference as `K`.  



{% highlight r %}
pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
{% endhighlight %}






and we use a harvest-based profit function with default parameters



{% highlight r %}
profit <- profit_harvest(price=1, c0 = 0.01) 
{% endhighlight %}




The `profit_harvest` function has the form \\( \Pi = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\).  Note that the R code defines a function from another function using a trick known as a _closure_.  Again we could write a custom profit function as long as it can take a vector stock size `x` and a scalar harvest level `h`.  Details for provided functions can be found in the manual, i.e. `?profit_harvest`. 


Now we must set up the discrete grids for stock size and havest levels (which will use same resolution as for stock), in order to calculate the SDP solution.   Here we set the gridsize to 100.  



{% highlight r %}
x_grid <- seq(0, 2 * K, length = 100)  
h_grid <- x_grid  
{% endhighlight %}





# Scenarios: 

We calculate the stochastic transition matrix for the probability of going from any state \\(x_t \\) to any other state \\(x_{t+1}\\) the following year, for each possible choice of harvest \\( h_t \\).  This provides a look-up table for the dynamic programming calculations.  


In the Sethi case, computing the distribution over multiple sources of noise is actually quite difficult.  Simulation turns out to be more efficient than numerically integrating over each distribution.  


## No Uncertainty 




{% highlight r %}
sigma_g <- 0.0    # Noise in population growth
sigma_m <- 0.0     # noise in stock assessment measurement
sigma_i <- 0.0     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
{% endhighlight %}



Find the transition matrix 



{% highlight r %}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
{% endhighlight %}




Find the optimum solution



{% highlight r %}
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
{% endhighlight %}




Simulate 



{% highlight r %}
sims_known <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
{% endhighlight %}




## Growth uncertainty 




{% highlight r %}
sigma_g <- 0.1    # Noise in population growth
sigma_m <- 0.0     # noise in stock assessment measurement
sigma_i <- 0.0     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
{% endhighlight %}



Find the transition matrix 



{% highlight r %}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
{% endhighlight %}




Find the optimum solution



{% highlight r %}
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
{% endhighlight %}




Simulate 



{% highlight r %}
sims_g <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
{% endhighlight %}





## Growth & stock measurement uncertainty 




{% highlight r %}
sigma_g <- 0.1    # Noise in population growth
sigma_m <- 0.1     # noise in stock assessment measurement
sigma_i <- 0.0     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
{% endhighlight %}




Find the transition matrix.  Use the simulation method to account for the extra uncertainties 



{% highlight r %}
require(snowfall) 
sfInit(parallel=TRUE, cpu=16)
{% endhighlight %}



{% highlight text %}
R Version:  R version 2.14.1 (2011-12-22) 

{% endhighlight %}



{% highlight r %}
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
{% endhighlight %}



{% highlight text %}
Library ggplot2 loaded.
{% endhighlight %}




Find the optimum solution



{% highlight r %}
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
{% endhighlight %}




Simulate 



{% highlight r %}
sims_gm <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
{% endhighlight %}







## Growth, stock measurement & implementation uncertainty 




{% highlight r %}
sigma_g <- 0.1    # Noise in population growth
sigma_m <- 0.1     # noise in stock assessment measurement
sigma_i <- 0.1     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
{% endhighlight %}




Find the transition matrix.  Use the simulation method to account for the extra uncertainties 



{% highlight r %}
require(snowfall) 
sfInit(parallel=TRUE, cpu=4)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
{% endhighlight %}



{% highlight text %}
Library ggplot2 loaded.
{% endhighlight %}




Find the optimum solution



{% highlight r %}
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
{% endhighlight %}




Simulate 



{% highlight r %}
sims_gmi <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
{% endhighlight %}






## Summarize and plot the results                                                   

R makes it easy to work with this big replicate data set.  We make data tidy (melt), fast (data.tables), and nicely labeled.



{% highlight r %}

sims <- list(known = sims_known, growth = sims_g, growth_stock = sims_gm, growth_stock_harvest = sims_gmi)

dat <- melt(sims, id=names(sims_known[[1]]))  
dt <- data.table(dat)
setnames(dt, c("L2", "L1"), c("reps", "uncertainty")) # names are nice
{% endhighlight %}




### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  



{% highlight r %}
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") + 
  facet_wrap(~uncertainty) 
{% endhighlight %}

![plot of chunk onerep](http://farm8.staticflickr.com/7235/7178234314_3dfe812f19_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.



{% highlight r %}
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2) + facet_wrap(~uncertainty)
{% endhighlight %}

![plot of chunk all](http://farm8.staticflickr.com/7082/7178235012_c3d3370758_o.png) 


We can also look at the harvest dynamics:



{% highlight r %}
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen") + facet_wrap(~uncertainty)
{% endhighlight %}

![plot of chunk harvestplot](http://farm8.staticflickr.com/7225/7178235688_676d514ea9_o.png) 


This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement: 



{% highlight r %}
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey") + facet_wrap(~uncertainty)
{% endhighlight %}

![plot of chunk escapement](http://farm9.staticflickr.com/8025/7178236184_b295b1ffed_o.png) 





{% highlight r %}
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, profit_fishing)) +
  geom_line(aes(time, policy_cost), col="darkblue") + facet_wrap(~uncertainty)
{% endhighlight %}



{% highlight text %}
Error: object 'profit_fishing' not found
{% endhighlight %}






{% highlight r %}
profits <-dt[ , sum(profit), by=c("reps", "uncertainty")] 
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
{% endhighlight %}

![plot of chunk unnamed-chunk-18](http://farm8.staticflickr.com/7092/7178236676_d0688d812b_o.png) 


Summary stats



{% highlight r %}
profits[, mean(V1), by=uncertainty]
{% endhighlight %}



{% highlight text %}
              uncertainty    V1
[1,]                known 10.03
[2,]               growth 33.48
[3,]         growth_stock 33.05
[4,] growth_stock_harvest 33.33
{% endhighlight %}



{% highlight r %}
profits[, sd(V1), by=uncertainty]
{% endhighlight %}



{% highlight text %}
              uncertainty    V1
[1,]                known 0.000
[2,]               growth 2.843
[3,]         growth_stock 3.036
[4,] growth_stock_harvest 2.827
{% endhighlight %}






# References



{% highlight text %}
Error: there is no package called 'knitcitations'
{% endhighlight %}



