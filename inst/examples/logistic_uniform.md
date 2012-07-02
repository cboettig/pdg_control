






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
c0 <- 0.01
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, c1 = c1)
```




with price = `1`, `c0` = `0.01` and `c1` = `0`. 




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







```r
scenario <- function(policy_g, policy_m, policy_i) {
    
    z_g <- function() 1 + (2 * runif(1, 0, 1) - 1) * policy_g
    z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * policy_m
    z_i <- function() 1 + (2 * runif(1, 0, 1) - 1) * policy_i
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps = 20000)
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

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, stock - x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7117/7479915294_00a0b3acea_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8012/7479915796_975bc84dd6_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7128/7479916322_492cc69e64_o.png) 




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

![plot of chunk onerep](http://farm9.staticflickr.com/8149/7479927898_681efb5f65_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7274/7479932860_c44a5d07f9_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8012/7479935244_1e4330e45b_o.png) 


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
<!-- Sun Jul  1 09:02:53 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 660.53 </TD> <TD align="right"> 674.15 </TD> <TD align="right"> 673.90 </TD> <TD align="right"> 674.15 </TD> <TD align="right"> 674.10 </TD> <TD align="right"> 674.05 </TD> <TD align="right"> 672.64 </TD> <TD align="right"> 673.82 </TD> <TD align="right"> 665.89 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 651.44 </TD> <TD align="right"> 666.29 </TD> <TD align="right"> 666.73 </TD> <TD align="right"> 668.55 </TD> <TD align="right"> 663.41 </TD> <TD align="right"> 662.26 </TD> <TD align="right"> 665.87 </TD> <TD align="right"> 668.17 </TD> <TD align="right"> 663.81 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 660.24 </TD> <TD align="right"> 684.29 </TD> <TD align="right"> 671.97 </TD> <TD align="right"> 686.14 </TD> <TD align="right"> 668.36 </TD> <TD align="right"> 664.38 </TD> <TD align="right"> 686.24 </TD> <TD align="right"> 672.70 </TD> <TD align="right"> 668.03 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 442.46 </TD> <TD align="right"> 561.65 </TD> <TD align="right"> 566.54 </TD> <TD align="right"> 585.23 </TD> <TD align="right"> 564.86 </TD> <TD align="right"> 587.01 </TD> <TD align="right"> 530.71 </TD> <TD align="right"> 584.27 </TD> <TD align="right"> 584.64 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 637.80 </TD> <TD align="right"> 652.41 </TD> <TD align="right"> 649.37 </TD> <TD align="right"> 652.63 </TD> <TD align="right"> 648.82 </TD> <TD align="right"> 650.09 </TD> <TD align="right"> 647.89 </TD> <TD align="right"> 650.26 </TD> <TD align="right"> 644.76 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 465.21 </TD> <TD align="right"> 559.37 </TD> <TD align="right"> 568.01 </TD> <TD align="right"> 585.51 </TD> <TD align="right"> 565.42 </TD> <TD align="right"> 584.46 </TD> <TD align="right"> 533.17 </TD> <TD align="right"> 592.40 </TD> <TD align="right"> 580.97 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 608.66 </TD> <TD align="right"> 647.40 </TD> <TD align="right"> 652.14 </TD> <TD align="right"> 637.49 </TD> <TD align="right"> 653.34 </TD> <TD align="right"> 638.67 </TD> <TD align="right"> 636.05 </TD> <TD align="right"> 631.30 </TD> <TD align="right"> 642.26 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 310.82 </TD> <TD align="right"> 406.54 </TD> <TD align="right"> 424.07 </TD> <TD align="right"> 483.51 </TD> <TD align="right"> 397.90 </TD> <TD align="right"> 462.25 </TD> <TD align="right"> 376.57 </TD> <TD align="right"> 473.13 </TD> <TD align="right"> 432.41 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 329.50 </TD> <TD align="right"> 405.50 </TD> <TD align="right"> 384.90 </TD> <TD align="right"> 440.31 </TD> <TD align="right"> 383.11 </TD> <TD align="right"> 468.04 </TD> <TD align="right"> 410.78 </TD> <TD align="right"> 426.65 </TD> <TD align="right"> 460.18 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Sun Jul  1 09:02:53 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 22.62 </TD> <TD align="right"> 20.37 </TD> <TD align="right"> 19.84 </TD> <TD align="right"> 19.57 </TD> <TD align="right"> 22.00 </TD> <TD align="right"> 22.03 </TD> <TD align="right"> 20.91 </TD> <TD align="right"> 21.14 </TD> <TD align="right"> 19.42 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 111.26 </TD> <TD align="right"> 120.16 </TD> <TD align="right"> 91.48 </TD> <TD align="right"> 109.47 </TD> <TD align="right"> 106.63 </TD> <TD align="right"> 103.53 </TD> <TD align="right"> 96.42 </TD> <TD align="right"> 102.63 </TD> <TD align="right"> 112.94 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 146.02 </TD> <TD align="right"> 58.08 </TD> <TD align="right"> 54.27 </TD> <TD align="right"> 20.37 </TD> <TD align="right"> 32.39 </TD> <TD align="right"> 24.57 </TD> <TD align="right"> 104.76 </TD> <TD align="right"> 19.44 </TD> <TD align="right"> 20.83 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 15.44 </TD> <TD align="right"> 11.74 </TD> <TD align="right"> 13.35 </TD> <TD align="right"> 11.70 </TD> <TD align="right"> 13.14 </TD> <TD align="right"> 12.34 </TD> <TD align="right"> 11.40 </TD> <TD align="right"> 12.45 </TD> <TD align="right"> 10.99 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 154.11 </TD> <TD align="right"> 121.86 </TD> <TD align="right"> 86.66 </TD> <TD align="right"> 93.88 </TD> <TD align="right"> 89.74 </TD> <TD align="right"> 90.55 </TD> <TD align="right"> 117.77 </TD> <TD align="right"> 97.08 </TD> <TD align="right"> 95.36 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 92.50 </TD> <TD align="right"> 93.88 </TD> <TD align="right"> 103.54 </TD> <TD align="right"> 94.77 </TD> <TD align="right"> 111.24 </TD> <TD align="right"> 95.34 </TD> <TD align="right"> 98.79 </TD> <TD align="right"> 94.10 </TD> <TD align="right"> 98.27 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 183.73 </TD> <TD align="right"> 180.20 </TD> <TD align="right"> 177.81 </TD> <TD align="right"> 142.32 </TD> <TD align="right"> 167.83 </TD> <TD align="right"> 161.60 </TD> <TD align="right"> 177.74 </TD> <TD align="right"> 151.09 </TD> <TD align="right"> 173.16 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 167.50 </TD> <TD align="right"> 182.49 </TD> <TD align="right"> 190.40 </TD> <TD align="right"> 176.95 </TD> <TD align="right"> 176.30 </TD> <TD align="right"> 176.42 </TD> <TD align="right"> 182.67 </TD> <TD align="right"> 169.48 </TD> <TD align="right"> 164.15 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


