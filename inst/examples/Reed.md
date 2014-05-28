`ro cache=FALSE, tidy=FALSE, warning=FALSE, comment=NA, message=FALSE, verbose=TRUE or`




# Reed Model

 * Author [Carl Boettiger](http://carlboettiger.info), <cboettig@gmail.com>
 * License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
 * Description:  Implements a numerical version of the SDP described in <a href="http://dx.doi.org/10.1016/0095-0696(79)90014-7">Reed (1979)</a>.


```
## Loading required package: pdgControl
```

```
## Loading required package: reshape2
```

```
## Loading required package: ggplot2
```

```
## Loading required package: data.table
```



Chose the state equation / population dynamics function


```r
f <- RickerAllee
K <- 10
pars <- c(2, K, 5)
```


We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  

<div> $$ \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x}, $$ </div> 


conditioned on h > x and x > 0,


```r
price <- 1
c0 <- 0
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, c1 = c1)
```


with price = 1, `c0` = 0 and `c1` = 0. 



```r
xmin <- 0
xmax <- K * 1.5
grid_n <- 200
```


We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from 0 to 15 on a grid of 200 points, and over an identical discrete grid of possible harvest values.  



```r
x_grid <- seq(xmin, xmax, length = grid_n)
h_grid <- x_grid
```




```r
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- 0.1
```


We will determine the optimal solution over a 25 time step window with boundary condition for stock at 0 and discounting rate of 0.05.  The Reed model considers a stochastic growth model 

<div> $$ x_{t+1} = z_g f(x_t) $$ </div> 

for the random variable `z_g`, given by 


```r
z_g <- function() 1 + rlnorm(1, 0, sigma_g)
pdfn <- function(P, s) dlnorm(P, 0, s)

# z_g <- function() 1+(2*runif(1, 0, 1)-1) * sigma_g pdfn <- function(P,
# s) dunif(P, 1 - s, 1 + s)
```






```r
p <- c(0.2, 7.5, 5)
SDP_Mat <- determine_SDP_matrix(f, p, x_grid, h_grid, sigma_g, pdfn)
```


### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward = 0)
opt$S
```

```
## [1] 7.161
```




Stationary optimal policy:  


```r
s_opt <- value_iteration(SDP_Mat, x_grid, h_grid, OptTime = 1000, xT, profit, 
    delta)
```



## Compare to Clark (noise-free)


```r
SDP_Mat2 <- determine_SDP_matrix(f, pars, x_grid, h_grid, 0.001, pdfn)
```


Bellman's algorithm to compute the optimal solution for all possible trajectories.


```r
det <- find_dp_optim(SDP_Mat2, x_grid, h_grid, OptTime, xT, profit, delta, reward = 0)
```




Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   


Plot the policy function (in terms of escapement, `x-h`, rather than harvest `h`) at equilibrium (first time-step):


```r
require(reshape2)
policies <- melt(data.frame(stock = x_grid, S = x_grid[opt$D[, 1]], D = x_grid[det$D[, 
    1]], Stationary = x_grid[s_opt$D]), id = "stock")
q1 <- ggplot(policies, aes(stock, stock - value, color = variable)) + geom_point(alpha = 0.5) + 
    xlab("stock size") + ylab("escapement")
q1
```

![plot of chunk policyfn_plot](http://farm4.staticflickr.com/3829/8856495917_0ee6c01d33_o.png) 


and the value function (at equilibrium):


```r
q2 <- qplot(x_grid, opt$V, xlab = "stock size", ylab = "value") + geom_vline(xintercept = opt$S)
q2
```

![plot of chunk valuefn_plot](http://farm4.staticflickr.com/3752/8857107876_d4a5ea9ed8_o.png) 






### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the Reed optimal harvest policy determined above.

No other sources of noise enter into the dynamics.  


```r
z_m <- function() 1
z_i <- function() 1
```




```r
sims <- lapply(1:100, function(i) {
    ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, opt$D, z_g, z_m, z_i)
})
```


The forward simulation algorithm needs an initial condition `x0` which we set equal to the carrying capacity, as well as our population dynamics `f`, parameters `pars`, grids, and noise coefficients.  Recall in the Reed case only `z_g`, growth, is stochastic.  


## Summarize and plot the results                                                   

R makes it easy to work with this big replicate data set.  We make data tidy (melt), fast (data.tables), and nicely labeled.


```r
dat <- melt(sims, id = names(sims[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "reps")  # names are nice
```


### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  


```r
p0 <- ggplot(subset(dt, reps == 1)) + geom_line(aes(time, fishstock)) + geom_abline(intercept = opt$S, 
    slope = 0) + geom_line(aes(time, harvest), col = "darkgreen")
p0
```

![plot of chunk p0](http://farm8.staticflickr.com/7429/8856499269_03790ea6f7_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 


```r
p1 <- ggplot(dt) + geom_abline(intercept = opt$S, slope = 0) + geom_abline(intercept = xT, 
    slope = 0, lty = 2)
p1 <- p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
p1
```

![plot of chunk p1](http://farm6.staticflickr.com/5330/8857111518_3f8f63dce1_o.png) 



# References






