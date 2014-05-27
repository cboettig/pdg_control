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
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps = 10000)
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, 
        reward = 0)
}
```




Determine the policies for each of the scenarios (noise combinations).



```r
lvl <- 0.5
```






```r
det <- scenario(0.01, 0, 0)
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

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, stock - x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8007/7688309276_7a2e980c53_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8166/7688309508_86156f1dcb_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8025/7688309746_d4f93e4cf4_o.png) 




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
noise <- list(det = c(0.01, 0, 0), low = c(0.1, 0.1, 0.1), growth = c(lvl, 
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

![plot of chunk onerep](http://farm9.staticflickr.com/8285/7688316766_e5d9a1cc52_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8009/7688319580_af03ef4b1c_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7116/7688321084_a48042c6c2_o.png) 


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
<!-- Tue Jul 31 19:01:07 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 675.11 </TD> <TD align="right"> 674.27 </TD> <TD align="right"> 673.02 </TD> <TD align="right"> 673.89 </TD> <TD align="right"> 674.23 </TD> <TD align="right"> 673.55 </TD> <TD align="right"> 671.76 </TD> <TD align="right"> 673.54 </TD> <TD align="right"> 671.73 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 665.83 </TD> <TD align="right"> 665.94 </TD> <TD align="right"> 663.64 </TD> <TD align="right"> 670.52 </TD> <TD align="right"> 667.69 </TD> <TD align="right"> 672.03 </TD> <TD align="right"> 664.73 </TD> <TD align="right"> 664.97 </TD> <TD align="right"> 667.88 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 689.76 </TD> <TD align="right"> 676.75 </TD> <TD align="right"> 663.14 </TD> <TD align="right"> 682.77 </TD> <TD align="right"> 689.30 </TD> <TD align="right"> 645.93 </TD> <TD align="right"> 668.38 </TD> <TD align="right"> 679.38 </TD> <TD align="right"> 669.94 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 560.21 </TD> <TD align="right"> 562.46 </TD> <TD align="right"> 565.81 </TD> <TD align="right"> 584.70 </TD> <TD align="right"> 565.79 </TD> <TD align="right"> 583.34 </TD> <TD align="right"> 545.29 </TD> <TD align="right"> 588.20 </TD> <TD align="right"> 587.20 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 650.93 </TD> <TD align="right"> 649.01 </TD> <TD align="right"> 649.86 </TD> <TD align="right"> 650.60 </TD> <TD align="right"> 648.50 </TD> <TD align="right"> 648.87 </TD> <TD align="right"> 648.67 </TD> <TD align="right"> 650.93 </TD> <TD align="right"> 649.83 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 549.23 </TD> <TD align="right"> 565.53 </TD> <TD align="right"> 560.99 </TD> <TD align="right"> 576.42 </TD> <TD align="right"> 560.15 </TD> <TD align="right"> 575.46 </TD> <TD align="right"> 511.63 </TD> <TD align="right"> 568.02 </TD> <TD align="right"> 560.16 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 642.03 </TD> <TD align="right"> 642.67 </TD> <TD align="right"> 656.75 </TD> <TD align="right"> 648.26 </TD> <TD align="right"> 628.22 </TD> <TD align="right"> 649.20 </TD> <TD align="right"> 632.45 </TD> <TD align="right"> 637.94 </TD> <TD align="right"> 640.23 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 368.54 </TD> <TD align="right"> 390.46 </TD> <TD align="right"> 402.25 </TD> <TD align="right"> 457.21 </TD> <TD align="right"> 446.00 </TD> <TD align="right"> 440.33 </TD> <TD align="right"> 393.53 </TD> <TD align="right"> 478.15 </TD> <TD align="right"> 450.59 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 395.65 </TD> <TD align="right"> 402.08 </TD> <TD align="right"> 415.58 </TD> <TD align="right"> 457.59 </TD> <TD align="right"> 381.31 </TD> <TD align="right"> 467.91 </TD> <TD align="right"> 370.44 </TD> <TD align="right"> 450.62 </TD> <TD align="right"> 434.50 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Tue Jul 31 19:01:07 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 1.96 </TD> <TD align="right"> 2.13 </TD> <TD align="right"> 2.11 </TD> <TD align="right"> 2.10 </TD> <TD align="right"> 2.01 </TD> <TD align="right"> 2.09 </TD> <TD align="right"> 1.93 </TD> <TD align="right"> 2.22 </TD> <TD align="right"> 2.21 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 20.63 </TD> <TD align="right"> 23.32 </TD> <TD align="right"> 20.01 </TD> <TD align="right"> 21.57 </TD> <TD align="right"> 20.67 </TD> <TD align="right"> 20.89 </TD> <TD align="right"> 19.90 </TD> <TD align="right"> 22.63 </TD> <TD align="right"> 22.31 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 102.67 </TD> <TD align="right"> 112.08 </TD> <TD align="right"> 107.00 </TD> <TD align="right"> 106.02 </TD> <TD align="right"> 100.44 </TD> <TD align="right"> 105.62 </TD> <TD align="right"> 99.65 </TD> <TD align="right"> 94.62 </TD> <TD align="right"> 99.94 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 55.88 </TD> <TD align="right"> 53.17 </TD> <TD align="right"> 33.63 </TD> <TD align="right"> 23.09 </TD> <TD align="right"> 56.45 </TD> <TD align="right"> 22.72 </TD> <TD align="right"> 85.41 </TD> <TD align="right"> 23.59 </TD> <TD align="right"> 18.10 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 13.02 </TD> <TD align="right"> 12.73 </TD> <TD align="right"> 12.33 </TD> <TD align="right"> 13.43 </TD> <TD align="right"> 12.23 </TD> <TD align="right"> 14.27 </TD> <TD align="right"> 11.91 </TD> <TD align="right"> 13.27 </TD> <TD align="right"> 12.12 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 142.79 </TD> <TD align="right"> 87.39 </TD> <TD align="right"> 89.12 </TD> <TD align="right"> 78.82 </TD> <TD align="right"> 120.40 </TD> <TD align="right"> 96.47 </TD> <TD align="right"> 152.08 </TD> <TD align="right"> 80.40 </TD> <TD align="right"> 88.43 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 95.16 </TD> <TD align="right"> 94.46 </TD> <TD align="right"> 88.30 </TD> <TD align="right"> 106.32 </TD> <TD align="right"> 101.46 </TD> <TD align="right"> 103.18 </TD> <TD align="right"> 100.09 </TD> <TD align="right"> 93.22 </TD> <TD align="right"> 92.74 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 179.96 </TD> <TD align="right"> 188.20 </TD> <TD align="right"> 177.03 </TD> <TD align="right"> 164.06 </TD> <TD align="right"> 163.62 </TD> <TD align="right"> 175.46 </TD> <TD align="right"> 178.39 </TD> <TD align="right"> 146.24 </TD> <TD align="right"> 161.24 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 179.56 </TD> <TD align="right"> 175.68 </TD> <TD align="right"> 179.91 </TD> <TD align="right"> 177.18 </TD> <TD align="right"> 175.81 </TD> <TD align="right"> 153.76 </TD> <TD align="right"> 181.18 </TD> <TD align="right"> 161.85 </TD> <TD align="right"> 164.32 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


