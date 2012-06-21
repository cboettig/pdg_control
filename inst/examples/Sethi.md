




# Sethi Model
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in 



Chose the state equation / population dynamics function as a logistic map: 



```r
f <- function(x,h,p){
	S = x - h
	p[1] * S * (1 - S/p[2]) + S
}
```




With parameters `r` = `1` and `K` = `100`.



```r
pars <- c(r, K)
```





We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  

<div> $$ \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x}, $$ </div> 

conditioned on h > x and x > 0,



```r
price <- 1
c0 <- 0.0
c1 <- 0
profit <- profit_harvest(price=price, c0 = c0, c1=c1) 
```




with price = `1`, `c0` = `0` and `c1` = `0`. 




```r
xmin <- 0
xmax <- 1.5 * K
grid_n <- 100
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `150` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




```r
x_grid <- seq(xmin, xmax, length = grid_n)  
h_grid <- x_grid  
```







```r
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- .5
sigma_m <- .5
sigma_i <- .5
```




We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  The Reed model considers a stochastic growth model 

<div> $$ x_{t+1} = z_g f(x_t) $$ </div> 

for the random variable `z_g`, given by 



```r
z_g <- function() 1+(2*runif(1, 0,  1)-1) * sigma_g
z_m <- function() 1+(2*runif(1, 0,  1)-1) * sigma_m
z_i <- function() 1+(2*runif(1, 0,  1)-1) * sigma_i
```




With `sigma_g` = `0.5`, `sigma_m` = `0.5`, `sigma_i` = `0.5`.


In the Sethi case, computing the distribution over multiple sources of noise is actually quite difficult.  Simulation turns out to be more efficient than numerically integrating over each distribution.  This code parallelizes the operation over four cores, but can be scaled to an arbitrary cluster. 



```r
require(snowfall) 
sfInit(parallel=TRUE, cpu=16)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=1e5)
```

```
Library ggplot2 loaded.
```



Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   



### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories. 



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```










```r
policy <- data.frame(stock = x_grid, value = opt$D[,1])
ggplot(policy) + 
  geom_point(aes(stock, stock-x_grid[value])) + 
	geom_smooth(aes(stock, stock-x_grid[value])) + ylab("escapement") 
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8153/7411017340_170939495e_o.png) 

```r

ggplot(policy) + 
  geom_point(aes(stock, x_grid[value])) + 
	geom_smooth(aes(stock, x_grid[value])) + ylab("harvest") 
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7110/7411017622_8f216dc141_o.png) 

```r

value <- data.frame(stock = x_grid, value=opt$V)
ggplot(value) + 
  geom_point(aes(stock, value)) +
  geom_smooth(aes(stock, value)) +
  ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7126/7411017892_9bf8f64ff5_o.png) 





### Simulate 

Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.



```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i)
})
```




The forward simulation algorithm needs an initial condition `x0` which we set equal to the carrying capacity, as well as our population dynamics `f`, parameters `pars`, grids, and noise coefficients.  Recall in the Reed case only `z_g`, growth, is stochastic, while in the Sethi example we now have three forms of stochasticity -- growth, measurement, and implementation.   


## Summarize and plot the results                                                   

R makes it easy to work with this big replicate data set.  We make data tidy (melt), fast (data.tables), and nicely labeled.



```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```





### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  



```r
p0 <- ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
p0
```

![plot of chunk p0](http://farm6.staticflickr.com/5112/7411018290_f0350e0f4d_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 



```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 <- p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
p1
```

![plot of chunk p1](http://farm6.staticflickr.com/5116/7411018680_e9fe9cc3d1_o.png) 





