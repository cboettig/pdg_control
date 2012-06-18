






# Calculating the value of information with highly nonlinear value function

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty, but under strong allee dynamics 




We consider the state dynamics 



```r
f <- RickerAllee
K <- 4
xT <- 2  # final value, also allee threshold
pars <- c(1.5, K, xT)
f
```

```
function (x, h, p) 
{
    sapply(x, function(x) {
        x <- max(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])
    })
}
<environment: namespace:pdgControl>
```




We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  \\( \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\),



```r
price <- 1
c0 <- 0.01
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, 
    c1 = c1)
```




with price `1`, `c0` `0.01` and `c1` `0`. 




```r
xmin <- 0
xmax <- 1.5 * K
grid_n <- 100
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `6` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




```r
x_grid <- seq(xmin, xmax, length = grid_n)
h_grid <- x_grid
```







```r
delta <- 0.05
xT <- 0
OptTime <- 25
```




We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  

# Scenarios: 

We use Monte Carlo integration over the noise processes to determine the transition matrix.  



```r
require(snowfall)
sfInit(parallel = TRUE, cpu = 16)
```

```
R Version:  R version 2.14.1 (2011-12-22) 

```







```r
scenario <- function(policy_g, policy_m, policy_i) {
    
    z_g <- function() rlnorm(1, 0, policy_g)
    z_m <- function() rlnorm(1, 0, policy_m)
    z_i <- function() rlnorm(1, 0, policy_i)
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, 
        z_g, z_m, z_i, reps = 20000)
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, 
        xT, profit, delta, reward = 0)
}
```






```r
det <- scenario(0, 0, 0)
```

```
Library ggplot2 loaded.
```

```r
g <- scenario(0.1, 0, 0)
```

```
Library ggplot2 loaded.
```

```r
m <- scenario(0, 0.1, 0)
```

```
Library ggplot2 loaded.
```

```r
i <- scenario(0, 0, 0.1)
```

```
Library ggplot2 loaded.
```

```r
gm <- scenario(0.1, 0.1, 0)
```

```
Library ggplot2 loaded.
```

```r
gi <- scenario(0.1, 0, 0.1)
```

```
Library ggplot2 loaded.
```

```r
gmi <- scenario(0.1, 0.1, 0.1)
```

```
Library ggplot2 loaded.
```






```r
det <- scenario(0.001, 0, 0)
```

```
Library ggplot2 loaded.
```




### plots



```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, det = det$D[, 
    1], g = g$D[, 1], m = m$D[, 1], gm = gm$D[, 1], gi = gi$D[, 
    1], gmi = gmi$D[, 1]), id = "stock")
ggplot(policy) + geom_point(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7093/7377022524_b8b775502d_o.png) 

```r
ggplot(policy) + geom_smooth(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7105/7377022796_713e9fd91c_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, gmi = gmi$V), 
    id = "stock")
ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8166/7377023094_6ae6a48f9f_o.png) 

```r
ggplot(value) + geom_smooth(aes(stock, value, 
    color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7222/7377023526_ffdd801eb5_o.png) 


## Simulations



```r
simulatereps <- function(opt, true_g, true_m, 
    true_i) {
    z_g <- function() rlnorm(1, 0, true_g)
    z_m <- function() rlnorm(1, 0, true_m)
    z_i <- function() rlnorm(1, 0, true_i)
    
    sims <- lapply(1:100, function(i) {
        ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, 
            opt$D, z_g, z_m, z_i, profit)
    })
    
    list(sims = sims, opt = opt, true_stochasticity = c(true_g, 
        true_m, true_i))
}
```






Corner cases 



```r
base <- simulatereps(det, 0, 0, 0)
reckless <- simulatereps(det, 0.1, 0.1, 0.1)
```






```r
full <- simulatereps(gmi, 0.1, 0.1, 0.1)
cautious <- simulatereps(gmi, 0, 0, 0)
```







```r
sims <- list(base = base$sims, full = full$sims, 
    reckless = reckless$sims, cautious = cautious$sims)
dat <- melt(sims, id = names(sims[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L2", "L1"), c("reps", "uncertainty"))  # names are nice
```




### Plots 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, 
    fishstock)) + geom_line(aes(time, harvest), col = "darkgreen") + 
    facet_wrap(~uncertainty)
```

![plot of chunk onerep](http://farm8.staticflickr.com/7072/7377024178_d189f694f3_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt)
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8145/7377024636_b57f130e4e_o.png) 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, 
    profit)) + facet_wrap(~uncertainty)
```

![The profits made in each time interval of a single replicate, by scenario](http://farm8.staticflickr.com/7083/7377025130_c0ceea9773_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7099/7377025530_0666ea70fe_o.png) 


Summary statistics 



```r
profits[, mean(V1), by = uncertainty]
```

```
     uncertainty     V1
[1,]        base 10.921
[2,]        full  9.731
[3,]    reckless  8.864
[4,]    cautious 10.676
```

```r
profits[, sd(V1), by = uncertainty]
```

```
     uncertainty    V1
[1,]        base 0.000
[2,]        full 1.763
[3,]    reckless 2.578
[4,]    cautious 0.000
```







# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


