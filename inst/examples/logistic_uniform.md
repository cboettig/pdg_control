






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

![plot of chunk sethiplots](http://farm9.staticflickr.com/8421/7498691408_b571722880_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8433/7498691854_181f66a860_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7126/7498692336_27021e0bb5_o.png) 




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

![plot of chunk onerep](http://farm9.staticflickr.com/8286/7498698468_0e71341759_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8020/7498701176_3f9df7dafd_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7128/7498702744_70f927dedb_o.png) 


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
<!-- Tue Jul  3 20:15:04 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 659.81 </TD> <TD align="right"> 673.34 </TD> <TD align="right"> 671.84 </TD> <TD align="right"> 673.35 </TD> <TD align="right"> 671.87 </TD> <TD align="right"> 673.36 </TD> <TD align="right"> 667.25 </TD> <TD align="right"> 671.84 </TD> <TD align="right"> 668.39 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 651.00 </TD> <TD align="right"> 669.18 </TD> <TD align="right"> 663.41 </TD> <TD align="right"> 669.85 </TD> <TD align="right"> 661.96 </TD> <TD align="right"> 663.92 </TD> <TD align="right"> 665.59 </TD> <TD align="right"> 660.94 </TD> <TD align="right"> 661.23 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 657.98 </TD> <TD align="right"> 679.43 </TD> <TD align="right"> 673.37 </TD> <TD align="right"> 689.91 </TD> <TD align="right"> 676.64 </TD> <TD align="right"> 654.06 </TD> <TD align="right"> 665.04 </TD> <TD align="right"> 663.19 </TD> <TD align="right"> 656.74 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 423.72 </TD> <TD align="right"> 558.53 </TD> <TD align="right"> 562.37 </TD> <TD align="right"> 584.57 </TD> <TD align="right"> 564.75 </TD> <TD align="right"> 586.09 </TD> <TD align="right"> 532.87 </TD> <TD align="right"> 584.08 </TD> <TD align="right"> 581.37 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 633.17 </TD> <TD align="right"> 649.15 </TD> <TD align="right"> 649.07 </TD> <TD align="right"> 649.89 </TD> <TD align="right"> 650.14 </TD> <TD align="right"> 648.05 </TD> <TD align="right"> 646.43 </TD> <TD align="right"> 648.64 </TD> <TD align="right"> 646.95 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 452.08 </TD> <TD align="right"> 557.95 </TD> <TD align="right"> 544.48 </TD> <TD align="right"> 571.24 </TD> <TD align="right"> 559.75 </TD> <TD align="right"> 584.32 </TD> <TD align="right"> 537.84 </TD> <TD align="right"> 568.27 </TD> <TD align="right"> 574.18 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 615.94 </TD> <TD align="right"> 631.32 </TD> <TD align="right"> 631.47 </TD> <TD align="right"> 644.91 </TD> <TD align="right"> 637.02 </TD> <TD align="right"> 626.62 </TD> <TD align="right"> 639.47 </TD> <TD align="right"> 638.83 </TD> <TD align="right"> 644.91 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 295.13 </TD> <TD align="right"> 376.18 </TD> <TD align="right"> 419.62 </TD> <TD align="right"> 457.61 </TD> <TD align="right"> 407.95 </TD> <TD align="right"> 470.97 </TD> <TD align="right"> 368.17 </TD> <TD align="right"> 474.43 </TD> <TD align="right"> 452.21 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 300.34 </TD> <TD align="right"> 390.99 </TD> <TD align="right"> 423.23 </TD> <TD align="right"> 448.29 </TD> <TD align="right"> 407.20 </TD> <TD align="right"> 426.56 </TD> <TD align="right"> 380.64 </TD> <TD align="right"> 467.14 </TD> <TD align="right"> 439.36 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Tue Jul  3 20:15:04 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 23.91 </TD> <TD align="right"> 20.34 </TD> <TD align="right"> 18.61 </TD> <TD align="right"> 20.42 </TD> <TD align="right"> 20.00 </TD> <TD align="right"> 19.92 </TD> <TD align="right"> 19.52 </TD> <TD align="right"> 21.47 </TD> <TD align="right"> 21.18 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 99.14 </TD> <TD align="right"> 120.78 </TD> <TD align="right"> 103.16 </TD> <TD align="right"> 102.11 </TD> <TD align="right"> 108.40 </TD> <TD align="right"> 93.39 </TD> <TD align="right"> 110.01 </TD> <TD align="right"> 91.05 </TD> <TD align="right"> 97.01 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 157.63 </TD> <TD align="right"> 74.68 </TD> <TD align="right"> 34.42 </TD> <TD align="right"> 20.78 </TD> <TD align="right"> 30.26 </TD> <TD align="right"> 23.99 </TD> <TD align="right"> 109.14 </TD> <TD align="right"> 24.04 </TD> <TD align="right"> 22.02 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 17.18 </TD> <TD align="right"> 13.16 </TD> <TD align="right"> 11.60 </TD> <TD align="right"> 13.46 </TD> <TD align="right"> 13.29 </TD> <TD align="right"> 13.67 </TD> <TD align="right"> 13.13 </TD> <TD align="right"> 12.20 </TD> <TD align="right"> 12.41 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 169.38 </TD> <TD align="right"> 91.12 </TD> <TD align="right"> 82.40 </TD> <TD align="right"> 92.57 </TD> <TD align="right"> 91.82 </TD> <TD align="right"> 103.30 </TD> <TD align="right"> 118.72 </TD> <TD align="right"> 89.16 </TD> <TD align="right"> 86.61 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 91.70 </TD> <TD align="right"> 85.86 </TD> <TD align="right"> 102.90 </TD> <TD align="right"> 106.21 </TD> <TD align="right"> 101.56 </TD> <TD align="right"> 104.08 </TD> <TD align="right"> 96.37 </TD> <TD align="right"> 101.63 </TD> <TD align="right"> 82.85 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 165.04 </TD> <TD align="right"> 171.77 </TD> <TD align="right"> 160.89 </TD> <TD align="right"> 166.68 </TD> <TD align="right"> 169.00 </TD> <TD align="right"> 147.65 </TD> <TD align="right"> 186.55 </TD> <TD align="right"> 165.02 </TD> <TD align="right"> 175.74 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 170.81 </TD> <TD align="right"> 176.10 </TD> <TD align="right"> 161.14 </TD> <TD align="right"> 156.53 </TD> <TD align="right"> 179.35 </TD> <TD align="right"> 182.02 </TD> <TD align="right"> 199.64 </TD> <TD align="right"> 148.99 </TD> <TD align="right"> 143.28 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


