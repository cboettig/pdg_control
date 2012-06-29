






# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   




Chose the state equation / population dynamics function



```r
f <- BevHolt
```





With parameters `A` = `1.5` and `B` = `0.05`.



```r
pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
```







We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  \\( \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\),



```r
price <- 1
c0 <- 0
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, 
    c1 = c1)
```




with price = `1`, `c0` = `0` and `c1` = `0`. 




```r
xmin <- 0
xmax <- 1.5 * K
grid_n <- 100
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `15` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




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
    
    z_g <- function() 1 + (2 * runif(1, 0, 1) - 1) * policy_g
    z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * policy_m
    z_i <- function() 1 + (2 * runif(1, 0, 1) - 1) * policy_i
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, 
        z_g, z_m, z_i, reps = 2e+05)
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, 
        xT, profit, delta, reward = 0)
}
```




Determine the policies for each of the scenarios (noise combinations).



```r
lvl <- 0.5
```






```r
det <- scenario(0, 0, 0)
```

```
Library ggplot2 loaded.
```






```r
all_low <- scenario(0.1, 0.1, 0.1)
```

```
Library ggplot2 loaded.
```






```r
g <- scenario(lvl, 0, 0)
```

```
Library ggplot2 loaded.
```






```r
m <- scenario(0, lvl, 0)
```

```
Library ggplot2 loaded.
```






```r
i <- scenario(0, 0, lvl)
```

```
Library ggplot2 loaded.
```






```r
gm <- scenario(lvl, lvl, 0)
```

```
Library ggplot2 loaded.
```






```r
gi <- scenario(lvl, 0, lvl)
```

```
Library ggplot2 loaded.
```






```r
mi <- scenario(0, lvl, lvl)
```

```
Library ggplot2 loaded.
```






```r
gmi <- scenario(lvl, lvl, lvl)
```

```
Library ggplot2 loaded.
```







```r
low <- all_low
```





### plots





```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, det = det$D[, 
    1], low = low$D[, 1], g = g$D[, 1], m = m$D[, 1], i = m$D[, 
    1], gm = gm$D[, 1], gi = gi$D[, 1], mi = mi$D[, 1], gmi = gmi$D[, 
    1]), id = "stock")

ggplot(policy) + geom_point(aes(stock, stock - 
    x_grid[value], color = variable)) + geom_smooth(aes(stock, 
    stock - x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8150/7466100344_a6c0c8208a_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], 
    color = variable)) + geom_smooth(aes(stock, x_grid[value], 
    color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8006/7466100652_10dd26c048_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    low = low$V, g = g$V, m = m$V, gm = gm$V, gi = gi$V, 
    mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7137/7466100932_9b8dcfba12_o.png) 


## Simulations



```r
simulatereps <- function(opt, true_g, true_m, 
    true_i) {
    
    z_g <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_g
    z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_m
    z_i <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_i
    
    sims <- lapply(1:100, function(i) {
        ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, 
            opt$D, z_g, z_m, z_i, profit)
    })
    
    sims
}
```





All cases



```r
policyfn <- list(det = det, low = low, g = g, 
    m = m, i = i, gm = gm, gi = gi, mi = mi, gmi = gmi)
noise <- list(det = c(0, 0, 0), low = c(0.1, 0.1, 
    0.1), growth = c(lvl, 0, 0), measure = c(0, lvl, 0), 
    implement = c(0, 0, lvl), growth_measure = c(lvl, lvl, 
        0), growth_implement = c(lvl, 0, lvl), measure_implement = c(0, 
        lvl, lvl), all = c(lvl, lvl, lvl))
allcases <- lapply(policyfn, function(policyfn_i) {
    lapply(noise, function(noise_i) {
        simulatereps(policyfn_i, noise_i[1], noise_i[2], 
            noise_i[3])
    })
})
```






```r
sims <- unlist(allcases, recursive = FALSE)
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

![plot of chunk onerep](http://farm8.staticflickr.com/7253/7466107252_c782545455_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7255/7466110580_99da161661_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7264/7466112530_a53c29ea95_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by = uncertainty]
sds <- profits[, sd(V1), by = uncertainty]
```






```r
require(xtable)
uncertainties <- names(noise)
print(xtable(matrix(means$V1, nrow = length(noise), 
    dimnames = list(uncertainties, uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Fri Jun 29 04:32:13 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 31.79 </TD> <TD align="right"> 33.18 </TD> <TD align="right"> 33.03 </TD> <TD align="right"> 32.09 </TD> <TD align="right"> 33.14 </TD> <TD align="right"> 32.07 </TD> <TD align="right"> 32.82 </TD> <TD align="right"> 31.97 </TD> <TD align="right"> 31.93 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 30.95 </TD> <TD align="right"> 33.03 </TD> <TD align="right"> 32.51 </TD> <TD align="right"> 32.04 </TD> <TD align="right"> 32.46 </TD> <TD align="right"> 32.15 </TD> <TD align="right"> 32.59 </TD> <TD align="right"> 31.94 </TD> <TD align="right"> 31.93 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 32.53 </TD> <TD align="right"> 32.39 </TD> <TD align="right"> 32.16 </TD> <TD align="right"> 31.93 </TD> <TD align="right"> 31.92 </TD> <TD align="right"> 32.30 </TD> <TD align="right"> 32.36 </TD> <TD align="right"> 32.27 </TD> <TD align="right"> 32.63 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 29.97 </TD> <TD align="right"> 29.43 </TD> <TD align="right"> 30.44 </TD> <TD align="right"> 30.77 </TD> <TD align="right"> 30.32 </TD> <TD align="right"> 30.94 </TD> <TD align="right"> 30.53 </TD> <TD align="right"> 30.84 </TD> <TD align="right"> 30.66 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 30.75 </TD> <TD align="right"> 32.31 </TD> <TD align="right"> 32.36 </TD> <TD align="right"> 32.06 </TD> <TD align="right"> 32.26 </TD> <TD align="right"> 32.06 </TD> <TD align="right"> 31.93 </TD> <TD align="right"> 31.98 </TD> <TD align="right"> 31.96 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 29.33 </TD> <TD align="right"> 27.98 </TD> <TD align="right"> 29.74 </TD> <TD align="right"> 27.84 </TD> <TD align="right"> 28.48 </TD> <TD align="right"> 29.03 </TD> <TD align="right"> 29.13 </TD> <TD align="right"> 30.29 </TD> <TD align="right"> 29.83 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 32.29 </TD> <TD align="right"> 30.24 </TD> <TD align="right"> 31.34 </TD> <TD align="right"> 31.38 </TD> <TD align="right"> 31.89 </TD> <TD align="right"> 32.45 </TD> <TD align="right"> 31.05 </TD> <TD align="right"> 32.04 </TD> <TD align="right"> 31.19 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 26.21 </TD> <TD align="right"> 25.77 </TD> <TD align="right"> 27.01 </TD> <TD align="right"> 28.33 </TD> <TD align="right"> 25.70 </TD> <TD align="right"> 28.01 </TD> <TD align="right"> 27.27 </TD> <TD align="right"> 29.02 </TD> <TD align="right"> 26.50 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 26.05 </TD> <TD align="right"> 23.25 </TD> <TD align="right"> 25.88 </TD> <TD align="right"> 27.21 </TD> <TD align="right"> 23.39 </TD> <TD align="right"> 26.26 </TD> <TD align="right"> 25.34 </TD> <TD align="right"> 26.52 </TD> <TD align="right"> 27.02 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), 
    dimnames = list(uncertainties, uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Fri Jun 29 04:32:13 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 1.93 </TD> <TD align="right"> 1.62 </TD> <TD align="right"> 1.59 </TD> <TD align="right"> 1.45 </TD> <TD align="right"> 1.50 </TD> <TD align="right"> 1.43 </TD> <TD align="right"> 1.81 </TD> <TD align="right"> 1.83 </TD> <TD align="right"> 1.59 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 8.09 </TD> <TD align="right"> 7.47 </TD> <TD align="right"> 7.33 </TD> <TD align="right"> 7.15 </TD> <TD align="right"> 7.11 </TD> <TD align="right"> 7.20 </TD> <TD align="right"> 7.99 </TD> <TD align="right"> 7.11 </TD> <TD align="right"> 7.56 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 2.30 </TD> <TD align="right"> 3.74 </TD> <TD align="right"> 1.45 </TD> <TD align="right"> 0.88 </TD> <TD align="right"> 2.41 </TD> <TD align="right"> 0.96 </TD> <TD align="right"> 1.32 </TD> <TD align="right"> 0.96 </TD> <TD align="right"> 0.97 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 1.01 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 0.79 </TD> <TD align="right"> 0.33 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 0.37 </TD> <TD align="right"> 0.97 </TD> <TD align="right"> 0.32 </TD> <TD align="right"> 0.40 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 8.05 </TD> <TD align="right"> 6.86 </TD> <TD align="right"> 7.27 </TD> <TD align="right"> 6.56 </TD> <TD align="right"> 7.35 </TD> <TD align="right"> 7.14 </TD> <TD align="right"> 6.67 </TD> <TD align="right"> 6.34 </TD> <TD align="right"> 7.48 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 7.87 </TD> <TD align="right"> 6.52 </TD> <TD align="right"> 7.73 </TD> <TD align="right"> 7.31 </TD> <TD align="right"> 8.66 </TD> <TD align="right"> 8.56 </TD> <TD align="right"> 7.61 </TD> <TD align="right"> 7.04 </TD> <TD align="right"> 7.52 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 7.17 </TD> <TD align="right"> 7.70 </TD> <TD align="right"> 6.63 </TD> <TD align="right"> 5.56 </TD> <TD align="right"> 7.63 </TD> <TD align="right"> 6.45 </TD> <TD align="right"> 6.86 </TD> <TD align="right"> 4.61 </TD> <TD align="right"> 7.05 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 7.93 </TD> <TD align="right"> 8.80 </TD> <TD align="right"> 8.73 </TD> <TD align="right"> 8.28 </TD> <TD align="right"> 8.12 </TD> <TD align="right"> 7.94 </TD> <TD align="right"> 8.21 </TD> <TD align="right"> 8.56 </TD> <TD align="right"> 8.17 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


