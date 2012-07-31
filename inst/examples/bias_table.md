






# Calculating the bias table 





Chose the state equation / population dynamics function



```r
f <- function(x, h, p) {
    sapply(x, function(x) {
        S = max(x - h, 0)
        p[1] * S * (1 - S/p[2]) + S
    })
}
```




With `K` = `100`, and variable choice of r.  


We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  \\( \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\),



```r
price <- 1
c0 <- 0
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, c1 = c1)
```




with price = `1`, `c0` = `0` and `c1` = `0`. 




```r
xmin <- 0
xmax <- 2.5 * K
grid_n <- 300
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `250` on a grid of `300` points, and over an identical discrete grid of possible harvest values.  




```r
x_grid <- seq(xmin, xmax, length = grid_n)
h_grid <- x_grid
```







```r
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- 0.01
sigma_m <- 0.5
```




We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  The Reed model considers a stochastic growth model 

<div> $$ x_{t+1} = z_g f(x_t) $$ </div> 

for the random variable `z_g`, given by 



```r
z_g <- function() 1 + (2 * runif(1, 0, 1) - 1) * sigma_g
z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * sigma_m
```




No other sources of noise enter into the dynamics.  



```r
z_i <- function() 1
```










## Scenario 1: low r

With parameter `r` = `1`



```r
pars1 <- c(r, K)
```






```r
require(snowfall)
sfInit(parallel = TRUE, cpu = 16)
```

```
R Version:  R version 2.14.1 (2011-12-22) 

```






```r
SDP_Mat1 <- SDP_by_simulation(f, pars1, x_grid, h_grid, z_g, z_m, 
    z_i, reps = 10000)
```

```
Library ggplot2 loaded.
```






```r
opt_low <- find_dp_optim(SDP_Mat1, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward = 0)
```





## Scenario 2: medium r

With parameter `r` = `1.5` 



```r
pars2 <- c(r, K)
```







```r
SDP_Mat2 <- SDP_by_simulation(f, pars2, x_grid, h_grid, z_g, z_m, 
    z_i, reps = 10000)
```

```
Library ggplot2 loaded.
```







```r
opt_med <- find_dp_optim(SDP_Mat2, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward = 0)
```




## Scenario 3: high r


With parameter `r` = `2` 



```r
pars3 <- c(r, K)
```






```r
SDP_Mat3 <- SDP_by_simulation(f, pars3, x_grid, h_grid, z_g, z_m, 
    z_i, reps = 10000)
```

```
Library ggplot2 loaded.
```






```r
opt_high <- find_dp_optim(SDP_Mat3, x_grid, h_grid, OptTime, xT, 
    profit, delta, reward = 0)
```





### plots



```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, low = opt_low$D[, 1], med = opt_med$D[, 
    1], high = opt_high$D[, 1]), id = "stock")

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, stock - x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7255/7681775570_fb4c0fe0c6_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8004/7681775794_392afa4afa_o.png) 

```r


value <- melt(data.frame(stock = x_grid, low = opt_low$V, med = opt_med$V, 
    high = opt_high$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8019/7681776006_58dfbfbb56_o.png) 




## Simulations



```r
simulatereps <- function(opt, pars) {
    sims <- lapply(1:5000, function(i) {
        ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, opt$D, z_g, z_m, z_i, 
            profit)
    })
    sims
}
```





All cases



```r
policyfn <- list(low = opt_low, med = opt_med, high = opt_high)
par_list <- list(pars1, pars2, pars3)

allcases <- lapply(policyfn, function(policyfn_i) {
    lapply(par_list, function(par) {
        simulatereps(policyfn_i, par)
    })
})
```






```r
sims <- unlist(allcases, recursive = FALSE)
dat <- melt(sims, id = names(sims[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L2", "L1"), c("reps", "parameter"))  # names are nice
```





### Plots 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, fishstock)) + 
    geom_line(aes(time, harvest), col = "darkgreen") + facet_wrap(~parameter)
```

![plot of chunk onerep](http://farm8.staticflickr.com/7107/7681821494_43678353b0_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~parameter)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8161/7681834216_bb3c7b188b_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "parameter")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~parameter)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7129/7681834788_651df17fe6_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by = parameter]
sds <- profits[, sd(V1), by = parameter]
```






```r
require(xtable)
scenarios <- c("low", "med", "high")
print(xtable(matrix(means$V1, nrow = length(scenarios), dimnames = list(scenarios, 
    scenarios))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Mon Jul 30 21:12:34 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> low </TH> <TH> med </TH> <TH> high </TH>  </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 583.71 </TD> <TD align="right"> 584.22 </TD> <TD align="right"> 584.36 </TD> </TR>
  <TR> <TD align="right"> med </TD> <TD align="right"> 787.21 </TD> <TD align="right"> 790.42 </TD> <TD align="right"> 792.23 </TD> </TR>
  <TR> <TD align="right"> high </TD> <TD align="right"> 965.00 </TD> <TD align="right"> 973.36 </TD> <TD align="right"> 976.40 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(scenarios), dimnames = list(scenarios, 
    scenarios))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Mon Jul 30 21:12:34 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> low </TH> <TH> med </TH> <TH> high </TH>  </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 23.33 </TD> <TD align="right"> 22.80 </TD> <TD align="right"> 21.93 </TD> </TR>
  <TR> <TD align="right"> med </TD> <TD align="right"> 39.38 </TD> <TD align="right"> 37.49 </TD> <TD align="right"> 36.88 </TD> </TR>
  <TR> <TD align="right"> high </TD> <TD align="right"> 60.76 </TD> <TD align="right"> 55.68 </TD> <TD align="right"> 54.52 </TD> </TR>
   </TABLE>





```r
k <- which.min(abs(K - x_grid))
data.frame(low = opt_low$V[k], med = opt_med$V[k], high = opt_high$V[k])
```

```
    low   med high
1 643.8 903.2 1142
```






