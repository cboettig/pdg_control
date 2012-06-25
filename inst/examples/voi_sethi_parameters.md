




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
profit <- profit_harvest(price = price, c0 = c0, 
    c1 = c1)
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
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, 
        z_g, z_m, z_i, reps = 20000)
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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7257/7441000118_07f9fd998a_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], 
    color = variable)) + geom_smooth(aes(stock, x_grid[value], 
    color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm6.staticflickr.com/5347/7441000648_0dfcd8c2f5_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    low = low$V, g = g$V, m = m$V, gm = gm$V, gi = gi$V, 
    mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7137/7441001286_f88dd0a87a_o.png) 


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

![plot of chunk onerep](http://farm9.staticflickr.com/8004/7441013866_77a8135163_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8163/7441020724_e380435bb5_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8151/7441024956_eaeb037482_o.png) 


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
<!-- Mon Jun 25 09:06:00 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> meas</TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 668.18 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 673.90 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 672.66 </TD> <TD align="right"> 672.73 </TD> <TD align="right"> 672.73 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 658.41 </TD> <TD align="right"> 667.17 </TD> <TD align="right"> 668.13 </TD> <TD align="right"> 668.70 </TD> <TD align="right"> 667.71 </TD> <TD align="right"> 667.11 </TD> <TD align="right"> 666.11 </TD> <TD align="right"> 665.49 </TD> <TD align="right"> 663.88 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 665.60 </TD> <TD align="right"> 658.17 </TD> <TD align="right"> 678.24 </TD> <TD align="right"> 672.08 </TD> <TD align="right"> 663.51 </TD> <TD align="right"> 663.35 </TD> <TD align="right"> 678.42 </TD> <TD align="right"> 673.02 </TD> <TD align="right"> 659.98 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 561.96 </TD> <TD align="right"> 563.93 </TD> <TD align="right"> 565.44 </TD> <TD align="right"> 580.35 </TD> <TD align="right"> 567.14 </TD> <TD align="right"> 586.47 </TD> <TD align="right"> 539.94 </TD> <TD align="right"> 585.76 </TD> <TD align="right"> 581.12 </TD> </TR>
  <TR> <TD align="right"> impl </TD> <TD align="right"> 642.66 </TD> <TD align="right"> 652.05 </TD> <TD align="right"> 648.74 </TD> <TD align="right"> 650.36 </TD> <TD align="right"> 648.10 </TD> <TD align="right"> 652.25 </TD> <TD align="right"> 648.68 </TD> <TD align="right"> 651.49 </TD> <TD align="right"> 647.17 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 570.61 </TD> <TD align="right"> 559.36 </TD> <TD align="right"> 549.87 </TD> <TD align="right"> 591.47 </TD> <TD align="right"> 546.31 </TD> <TD align="right"> 568.15 </TD> <TD align="right"> 526.44 </TD> <TD align="right"> 571.87 </TD> <TD align="right"> 570.54 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 656.30 </TD> <TD align="right"> 630.27 </TD> <TD align="right"> 631.36 </TD> <TD align="right"> 649.04 </TD> <TD align="right"> 642.05 </TD> <TD align="right"> 636.24 </TD> <TD align="right"> 626.19 </TD> <TD align="right"> 644.27 </TD> <TD align="right"> 629.62 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 409.87 </TD> <TD align="right"> 414.02 </TD> <TD align="right"> 409.90 </TD> <TD align="right"> 453.13 </TD> <TD align="right"> 390.19 </TD> <TD align="right"> 454.33 </TD> <TD align="right"> 383.75 </TD> <TD align="right"> 461.56 </TD> <TD align="right"> 437.81 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 414.06 </TD> <TD align="right"> 414.96 </TD> <TD align="right"> 380.41 </TD> <TD align="right"> 444.25 </TD> <TD align="right"> 407.72 </TD> <TD align="right"> 471.61 </TD> <TD align="right"> 374.73 </TD> <TD align="right"> 463.12 </TD> <TD align="right"> 449.96 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), 
    dimnames = list(uncertainties, uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Mon Jun 25 09:06:00 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> meas </TH> <TH> impl </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 22.80 </TD> <TD align="right"> 19.07 </TD> <TD align="right"> 22.03 </TD> <TD align="right"> 17.68 </TD> <TD align="right"> 19.26 </TD> <TD align="right"> 22.90 </TD> <TD align="right"> 17.36 </TD> <TD align="right"> 21.15 </TD> <TD align="right"> 22.38 </TD> </TR>
  <TR> <TD align="right"> grow </TD> <TD align="right"> 106.87 </TD> <TD align="right"> 107.47 </TD> <TD align="right"> 100.07 </TD> <TD align="right"> 110.41 </TD> <TD align="right"> 107.83 </TD> <TD align="right"> 93.08 </TD> <TD align="right"> 110.18 </TD> <TD align="right"> 119.66 </TD> <TD align="right"> 98.79 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 55.92 </TD> <TD align="right"> 61.69 </TD> <TD align="right"> 31.18 </TD> <TD align="right"> 25.07 </TD> <TD align="right"> 29.78 </TD> <TD align="right"> 20.51 </TD> <TD align="right"> 102.15 </TD> <TD align="right"> 24.69 </TD> <TD align="right"> 18.85 </TD> </TR>
  <TR> <TD align="right"> impl </TD> <TD align="right"> 12.52 </TD> <TD align="right"> 11.34 </TD> <TD align="right"> 13.38 </TD> <TD align="right"> 12.56 </TD> <TD align="right"> 13.01 </TD> <TD align="right"> 13.27 </TD> <TD align="right"> 11.40 </TD> <TD align="right"> 14.43 </TD> <TD align="right"> 11.62 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 106.36 </TD> <TD align="right"> 96.71 </TD> <TD align="right"> 92.00 </TD> <TD align="right"> 94.33 </TD> <TD align="right"> 127.93 </TD> <TD align="right"> 92.74 </TD> <TD align="right"> 137.62 </TD> <TD align="right"> 99.40 </TD> <TD align="right"> 82.45 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 93.77 </TD> <TD align="right"> 102.25 </TD> <TD align="right"> 106.51 </TD> <TD align="right"> 96.90 </TD> <TD align="right"> 102.85 </TD> <TD align="right"> 97.06 </TD> <TD align="right"> 79.78 </TD> <TD align="right"> 98.75 </TD> <TD align="right"> 97.92 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 178.13 </TD> <TD align="right"> 168.77 </TD> <TD align="right"> 165.51 </TD> <TD align="right"> 162.55 </TD> <TD align="right"> 179.85 </TD> <TD align="right"> 167.82 </TD> <TD align="right"> 175.59 </TD> <TD align="right"> 154.52 </TD> <TD align="right"> 155.50 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 166.82 </TD> <TD align="right"> 179.12 </TD> <TD align="right"> 189.97 </TD> <TD align="right"> 191.33 </TD> <TD align="right"> 195.58 </TD> <TD align="right"> 169.33 </TD> <TD align="right"> 182.70 </TD> <TD align="right"> 163.87 </TD> <TD align="right"> 162.37 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


