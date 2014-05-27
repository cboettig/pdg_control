---
layout: page
---








# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   




Chose the state equation / population dynamics function



```r
f <- function(x, h, p) {
    sapply(x, function(x) {
        S = max(x - h, 0)
        p[1] * S * (1 - S/p[2]) + S
    })
}
```




With parameters `r` = `1` and `K` = `100`.



```r
pars <- c(r, K)
```




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
xmax <- 1.5 * K
grid_n <- 500
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `150` on a grid of `500` points, and over an identical discrete grid of possible harvest values.  




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
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps = 1000)
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, 
        reward = 0)
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
policy <- melt(data.frame(stock = x_grid, det = det$D[, 1], low = low$D[, 
    1], g = g$D[, 1], m = m$D[, 1], i = m$D[, 1], gm = gm$D[, 1], gi = gi$D[, 
    1], mi = mi$D[, 1], gmi = gmi$D[, 1]), id = "stock")

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, stock - x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8008/7465381326_82782c0832_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, x_grid[value], color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8008/7465381782_ca1b81ec65_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8164/7465382200_4cd4be7aea_o.png) 


## Simulations



```r
simulatereps <- function(opt, true_g, true_m, true_i) {
    
    z_g <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_g
    z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_m
    z_i <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_i
    
    sims <- lapply(1:100, function(i) {
        ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, opt$D, z_g, z_m, z_i, 
            profit)
    })
    
    sims
}
```





All cases



```r
policyfn <- list(det = det, low = low, g = g, m = m, i = i, gm = gm, 
    gi = gi, mi = mi, gmi = gmi)
noise <- list(det = c(0, 0, 0), low = c(0.1, 0.1, 0.1), growth = c(lvl, 
    0, 0), measure = c(0, lvl, 0), implement = c(0, 0, lvl), growth_measure = c(lvl, 
    lvl, 0), growth_implement = c(lvl, 0, lvl), measure_implement = c(0, lvl, 
    lvl), all = c(lvl, lvl, lvl))
allcases <- lapply(policyfn, function(policyfn_i) {
    lapply(noise, function(noise_i) {
        simulatereps(policyfn_i, noise_i[1], noise_i[2], noise_i[3])
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
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, fishstock)) + 
    geom_line(aes(time, harvest), col = "darkgreen") + facet_wrap(~uncertainty)
```

![plot of chunk onerep](http://farm9.staticflickr.com/8025/7465424810_85b6ff7fbf_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7246/7465440122_a2ce5e4523_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8007/7465448150_44a724c88a_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by = uncertainty]
sds <- profits[, sd(V1), by = uncertainty]
```






```r
require(xtable)
uncertainties <- names(noise)
print(xtable(matrix(means$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Fri Jun 29 00:50:34 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 674.55 </TD> <TD align="right"> 673.65 </TD> <TD align="right"> 670.79 </TD> <TD align="right"> 674.55 </TD> <TD align="right"> 672.44 </TD> <TD align="right"> 673.57 </TD> <TD align="right"> 672.14 </TD> <TD align="right"> 672.75 </TD> <TD align="right"> 670.64 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 666.93 </TD> <TD align="right"> 665.62 </TD> <TD align="right"> 663.05 </TD> <TD align="right"> 668.57 </TD> <TD align="right"> 662.47 </TD> <TD align="right"> 668.70 </TD> <TD align="right"> 658.65 </TD> <TD align="right"> 663.27 </TD> <TD align="right"> 661.07 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 668.38 </TD> <TD align="right"> 678.18 </TD> <TD align="right"> 657.48 </TD> <TD align="right"> 648.18 </TD> <TD align="right"> 661.12 </TD> <TD align="right"> 679.13 </TD> <TD align="right"> 675.45 </TD> <TD align="right"> 670.19 </TD> <TD align="right"> 643.55 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 572.00 </TD> <TD align="right"> 560.38 </TD> <TD align="right"> 563.81 </TD> <TD align="right"> 580.61 </TD> <TD align="right"> 538.48 </TD> <TD align="right"> 582.72 </TD> <TD align="right"> 535.41 </TD> <TD align="right"> 581.51 </TD> <TD align="right"> 582.02 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 650.25 </TD> <TD align="right"> 650.09 </TD> <TD align="right"> 647.11 </TD> <TD align="right"> 650.86 </TD> <TD align="right"> 648.54 </TD> <TD align="right"> 650.84 </TD> <TD align="right"> 646.00 </TD> <TD align="right"> 647.16 </TD> <TD align="right"> 645.22 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 562.86 </TD> <TD align="right"> 546.40 </TD> <TD align="right"> 554.54 </TD> <TD align="right"> 569.10 </TD> <TD align="right"> 568.16 </TD> <TD align="right"> 570.31 </TD> <TD align="right"> 500.90 </TD> <TD align="right"> 569.37 </TD> <TD align="right"> 579.61 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 655.75 </TD> <TD align="right"> 651.53 </TD> <TD align="right"> 636.84 </TD> <TD align="right"> 630.30 </TD> <TD align="right"> 636.91 </TD> <TD align="right"> 633.56 </TD> <TD align="right"> 627.58 </TD> <TD align="right"> 639.51 </TD> <TD align="right"> 637.44 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 404.38 </TD> <TD align="right"> 388.77 </TD> <TD align="right"> 425.63 </TD> <TD align="right"> 395.62 </TD> <TD align="right"> 415.00 </TD> <TD align="right"> 483.36 </TD> <TD align="right"> 372.95 </TD> <TD align="right"> 454.57 </TD> <TD align="right"> 399.73 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 424.08 </TD> <TD align="right"> 419.80 </TD> <TD align="right"> 394.17 </TD> <TD align="right"> 440.61 </TD> <TD align="right"> 431.73 </TD> <TD align="right"> 465.14 </TD> <TD align="right"> 398.18 </TD> <TD align="right"> 458.54 </TD> <TD align="right"> 448.91 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Fri Jun 29 00:50:34 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 22.24 </TD> <TD align="right"> 20.61 </TD> <TD align="right"> 23.36 </TD> <TD align="right"> 21.64 </TD> <TD align="right"> 22.58 </TD> <TD align="right"> 24.35 </TD> <TD align="right"> 15.16 </TD> <TD align="right"> 22.05 </TD> <TD align="right"> 21.87 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 96.70 </TD> <TD align="right"> 97.71 </TD> <TD align="right"> 114.29 </TD> <TD align="right"> 102.16 </TD> <TD align="right"> 109.00 </TD> <TD align="right"> 92.29 </TD> <TD align="right"> 90.72 </TD> <TD align="right"> 111.55 </TD> <TD align="right"> 90.67 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 25.16 </TD> <TD align="right"> 60.40 </TD> <TD align="right"> 47.15 </TD> <TD align="right"> 28.83 </TD> <TD align="right"> 106.92 </TD> <TD align="right"> 25.19 </TD> <TD align="right"> 103.13 </TD> <TD align="right"> 23.15 </TD> <TD align="right"> 21.56 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 12.79 </TD> <TD align="right"> 12.85 </TD> <TD align="right"> 13.06 </TD> <TD align="right"> 12.51 </TD> <TD align="right"> 13.18 </TD> <TD align="right"> 12.02 </TD> <TD align="right"> 12.37 </TD> <TD align="right"> 14.22 </TD> <TD align="right"> 12.43 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 92.26 </TD> <TD align="right"> 102.60 </TD> <TD align="right"> 82.59 </TD> <TD align="right"> 95.05 </TD> <TD align="right"> 96.76 </TD> <TD align="right"> 81.68 </TD> <TD align="right"> 141.29 </TD> <TD align="right"> 93.09 </TD> <TD align="right"> 82.26 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 85.44 </TD> <TD align="right"> 100.55 </TD> <TD align="right"> 93.55 </TD> <TD align="right"> 89.28 </TD> <TD align="right"> 104.50 </TD> <TD align="right"> 105.57 </TD> <TD align="right"> 92.68 </TD> <TD align="right"> 110.01 </TD> <TD align="right"> 99.00 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 180.95 </TD> <TD align="right"> 190.16 </TD> <TD align="right"> 164.35 </TD> <TD align="right"> 182.05 </TD> <TD align="right"> 166.10 </TD> <TD align="right"> 149.17 </TD> <TD align="right"> 186.60 </TD> <TD align="right"> 159.40 </TD> <TD align="right"> 173.12 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 166.84 </TD> <TD align="right"> 179.09 </TD> <TD align="right"> 181.99 </TD> <TD align="right"> 188.21 </TD> <TD align="right"> 182.37 </TD> <TD align="right"> 177.96 </TD> <TD align="right"> 172.24 </TD> <TD align="right"> 151.33 </TD> <TD align="right"> 173.53 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


