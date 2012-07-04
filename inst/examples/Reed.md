




# Reed Model

 * Author [Carl Boettiger](http://carlboettiger.info), <cboettig@gmail.com>
 * License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
 * Description:  Implements a numerical version of the SDP described in Reed, (1979).





Chose the state equation / population dynamics function



```r
f <- BevHolt 
```




With parameters `A` = `1.5` and `B` = `0.005`.



```r
pars <- c(A, B)
```




We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  

<div> $$ \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x}, $$ </div> 

conditioned on h > x and x > 0,



```r
price <- 1
c0 <- 5
c1 <- 0
profit <- profit_harvest(price=price, c0 = c0, c1=c1) 
```




with price = `1`, `c0` = `5` and `c1` = `0`. 




```r
xmin <- 0
xmax <- 1.5 * K
```

```
Error: object 'K' not found
```

```r
grid_n <- 100
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to 

```

Error in eval(expr, envir, enclos) : object 'xmax' not found

```

 on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




```r
x_grid <- seq(xmin, xmax, length = grid_n)  
```

```
Error: object 'xmax' not found
```

```r
h_grid <- x_grid  
```

```
Error: object 'x_grid' not found
```







```r
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- .2
```




We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  The Reed model considers a stochastic growth model 

<div> $$ x_{t+1} = z_g f(x_t) $$ </div> 

for the random variable `z_g`, given by 



```r
z_g <- function() rlnorm(1,0, sigma_g)  # 1+(2*runif(1, 0,  1)-1) * sigma_g
```




No other sources of noise enter into the dynamics.  



```r
z_m <- function() 1
z_i <- function() 1
```








```r
pdfn <- function(P, s){
  # dunif(P, 1 - s, 1 + s)
  dlnorm(P, 0, s)
}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g, pdfn)
```

```
Error: object 'x_grid' not found
```





### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories.



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward=0)
```

```
Error: object 'x_grid' not found
```





## Compare to Clark (noise-free)



```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, 0.01, pdfn2)
```

```
Error: object 'pdfn2' not found
```




Bellman's algorithm to compute the optimal solution for all possible trajectories.



```r
det <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward=0)
```

```
Error: object 'x_grid' not found
```






Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   


Plot the policy function (in terms of escapement, `x-h`, rather than harvest `h`) at equilibrium (first time-step):



```r
require(reshape2)
policies <- melt(data.frame(stock=x_grid, S = x_grid[opt$D[,1]], D = x_grid[det$D[,1]]), id="stock")
```

```
Error: object 'x_grid' not found
```

```r
q1 <- ggplot(policies, aes(stock, stock - value, color=variable)) + geom_point() + xlab("stock size") + ylab("escapement") 
```

```
Error: object 'policies' not found
```

```r
q1
```

```
Error: object 'q1' not found
```




and the value function (at equilibrium):



```r
q2 <- qplot(x_grid, opt$V, xlab="stock size", ylab="value") + 
geom_vline(xintercept=opt$S)
```

```
Error: object 'opt' not found
```

```r
q2
```

```
Error: object 'q2' not found
```








### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the Reed optimal harvest policy determined above.



```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i)
})
```

```
Error: object 'opt' not found
```




The forward simulation algorithm needs an initial condition `x0` which we set equal to the carrying capacity, as well as our population dynamics `f`, parameters `pars`, grids, and noise coefficients.  Recall in the Reed case only `z_g`, growth, is stochastic.  


## Summarize and plot the results                                                   

R makes it easy to work with this big replicate data set.  We make data tidy (melt), fast (data.tables), and nicely labeled.



```r
dat <- melt(sims, id=names(sims[[1]]))  
```

```
Error: object 'sims' not found
```

```r
dt <- data.table(dat)
```

```
Error: object 'dat' not found
```

```r
setnames(dt, "L1", "reps") # names are nice
```

```
Error: x is not a data.table
```




### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  



```r
p0 <- ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
```

```
Error: object 'reps' not found
```

```r
p0
```

```
Error: object 'p0' not found
```





This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 



```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
```

```
Error: ggplot2 doesn't know how to deal with data of class function
```

```r
p1 <- p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

```
Error: object 'p1' not found
```

```r
p1
```

```
Error: object 'p1' not found
```





# References

<p>Reed WJ (1979).
&ldquo;Optimal Escapement Levels in Stochastic And Deterministic Harvesting Models.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>6</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/0095-0696(79)90014-7">http://dx.doi.org/10.1016/0095-0696(79)90014-7</a>.






```r
options(device=orig)
```

```
Error: object 'orig' not found
```



