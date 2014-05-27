---
layout: page
---








# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   




Chose the state equation / population dynamics function



```r
f <- BevHolt
```





With parameters `A` = `1.5` and `B` = `0.005`.



```r
pars <- c(A, B)
K <- (pars[1] - 1)/pars[2]
```





We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  \\( \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\),



```r
price <- 1
c0 <- 0.1
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, c1 = c1)
```




with price = `1`, `c0` = `0.1` and `c1` = `0`. 




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
    
    z_g <- function() rlnorm(1, 0, policy_g)
    z_m <- function() rlnorm(1, 0, policy_m)
    z_i <- function() rlnorm(1, 0, policy_i)
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps = 10000)
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, 
        reward = 0)
}
```




Determine the policies for each of the scenarios (noise combinations).



```r
lvl <- 0.2
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

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, stock - x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8002/7496376204_bd4e755a5d_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7276/7496376608_646eb972ee_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7119/7496377030_9c116a7121_o.png) 




## Simulations



```r
simulatereps <- function(opt, true_g, true_m, true_i) {
    
    z_g <- function() rlnorm(1, 0, true_g)
    z_m <- function() rlnorm(1, 0, true_m)
    z_i <- function() rlnorm(1, 0, true_i)
    
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
noise <- list(det = c(0, 0, 0), low = c(0.1, 0.1, 0.1), gro = c(lvl, 
    0, 0), meas = c(0, lvl, 0), imp = c(0, 0, lvl), gro_meas = c(lvl, lvl, 0), 
    gro_imp = c(lvl, 0, lvl), meas_imp = c(0, lvl, lvl), all = c(lvl, lvl, lvl))
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

![plot of chunk onerep](http://farm8.staticflickr.com/7252/7496385592_773f59e177_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8291/7496389552_9f152edcfc_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8283/7496391508_fd2133e2af_o.png) 


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
<!-- Tue Jul  3 12:36:01 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 316.24 </TD> <TD align="right"> 330.66 </TD> <TD align="right"> 329.06 </TD> <TD align="right"> 331.27 </TD> <TD align="right"> 329.77 </TD> <TD align="right"> 329.77 </TD> <TD align="right"> 326.77 </TD> <TD align="right"> 329.57 </TD> <TD align="right"> 326.78 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 301.61 </TD> <TD align="right"> 329.57 </TD> <TD align="right"> 327.78 </TD> <TD align="right"> 331.68 </TD> <TD align="right"> 333.45 </TD> <TD align="right"> 331.77 </TD> <TD align="right"> 334.72 </TD> <TD align="right"> 333.35 </TD> <TD align="right"> 327.86 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 348.19 </TD> <TD align="right"> 352.67 </TD> <TD align="right"> 354.21 </TD> <TD align="right"> 357.24 </TD> <TD align="right"> 357.61 </TD> <TD align="right"> 355.04 </TD> <TD align="right"> 367.41 </TD> <TD align="right"> 358.06 </TD> <TD align="right"> 367.83 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 253.33 </TD> <TD align="right"> 313.06 </TD> <TD align="right"> 298.50 </TD> <TD align="right"> 315.11 </TD> <TD align="right"> 304.88 </TD> <TD align="right"> 317.00 </TD> <TD align="right"> 309.71 </TD> <TD align="right"> 313.31 </TD> <TD align="right"> 315.40 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 302.83 </TD> <TD align="right"> 325.78 </TD> <TD align="right"> 323.23 </TD> <TD align="right"> 325.82 </TD> <TD align="right"> 323.10 </TD> <TD align="right"> 324.92 </TD> <TD align="right"> 323.46 </TD> <TD align="right"> 323.98 </TD> <TD align="right"> 322.16 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 257.21 </TD> <TD align="right"> 337.69 </TD> <TD align="right"> 323.83 </TD> <TD align="right"> 336.66 </TD> <TD align="right"> 333.21 </TD> <TD align="right"> 335.82 </TD> <TD align="right"> 339.23 </TD> <TD align="right"> 335.03 </TD> <TD align="right"> 341.52 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 310.30 </TD> <TD align="right"> 353.87 </TD> <TD align="right"> 351.70 </TD> <TD align="right"> 352.72 </TD> <TD align="right"> 345.56 </TD> <TD align="right"> 352.21 </TD> <TD align="right"> 334.11 </TD> <TD align="right"> 352.50 </TD> <TD align="right"> 355.34 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 217.61 </TD> <TD align="right"> 296.58 </TD> <TD align="right"> 288.99 </TD> <TD align="right"> 301.51 </TD> <TD align="right"> 283.85 </TD> <TD align="right"> 304.61 </TD> <TD align="right"> 294.43 </TD> <TD align="right"> 303.46 </TD> <TD align="right"> 305.55 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 240.35 </TD> <TD align="right"> 309.04 </TD> <TD align="right"> 296.89 </TD> <TD align="right"> 327.78 </TD> <TD align="right"> 314.79 </TD> <TD align="right"> 318.22 </TD> <TD align="right"> 322.38 </TD> <TD align="right"> 332.91 </TD> <TD align="right"> 323.74 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Tue Jul  3 12:36:01 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 57.37 </TD> <TD align="right"> 34.82 </TD> <TD align="right"> 25.37 </TD> <TD align="right"> 27.15 </TD> <TD align="right"> 26.23 </TD> <TD align="right"> 28.31 </TD> <TD align="right"> 30.56 </TD> <TD align="right"> 26.79 </TD> <TD align="right"> 29.64 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 50.96 </TD> <TD align="right"> 56.84 </TD> <TD align="right"> 66.58 </TD> <TD align="right"> 51.05 </TD> <TD align="right"> 61.20 </TD> <TD align="right"> 62.30 </TD> <TD align="right"> 57.90 </TD> <TD align="right"> 60.30 </TD> <TD align="right"> 62.76 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 68.82 </TD> <TD align="right"> 12.81 </TD> <TD align="right"> 48.18 </TD> <TD align="right"> 9.69 </TD> <TD align="right"> 39.74 </TD> <TD align="right"> 8.46 </TD> <TD align="right"> 30.40 </TD> <TD align="right"> 23.51 </TD> <TD align="right"> 8.89 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 37.56 </TD> <TD align="right"> 5.65 </TD> <TD align="right"> 6.21 </TD> <TD align="right"> 5.83 </TD> <TD align="right"> 6.93 </TD> <TD align="right"> 7.04 </TD> <TD align="right"> 5.56 </TD> <TD align="right"> 5.95 </TD> <TD align="right"> 6.59 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 83.30 </TD> <TD align="right"> 56.98 </TD> <TD align="right"> 55.09 </TD> <TD align="right"> 51.29 </TD> <TD align="right"> 66.41 </TD> <TD align="right"> 50.17 </TD> <TD align="right"> 66.77 </TD> <TD align="right"> 53.32 </TD> <TD align="right"> 46.86 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 71.95 </TD> <TD align="right"> 60.21 </TD> <TD align="right"> 57.63 </TD> <TD align="right"> 53.31 </TD> <TD align="right"> 54.94 </TD> <TD align="right"> 62.93 </TD> <TD align="right"> 62.52 </TD> <TD align="right"> 59.57 </TD> <TD align="right"> 59.50 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 74.11 </TD> <TD align="right"> 44.48 </TD> <TD align="right"> 61.14 </TD> <TD align="right"> 43.30 </TD> <TD align="right"> 66.34 </TD> <TD align="right"> 39.62 </TD> <TD align="right"> 48.02 </TD> <TD align="right"> 38.42 </TD> <TD align="right"> 32.77 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 80.88 </TD> <TD align="right"> 82.39 </TD> <TD align="right"> 77.43 </TD> <TD align="right"> 60.41 </TD> <TD align="right"> 63.46 </TD> <TD align="right"> 68.65 </TD> <TD align="right"> 72.09 </TD> <TD align="right"> 61.63 </TD> <TD align="right"> 70.74 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


