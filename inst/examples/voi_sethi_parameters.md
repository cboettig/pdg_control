






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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7137/7456298040_1f4ea25cb8_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], 
    color = variable)) + geom_smooth(aes(stock, x_grid[value], 
    color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7270/7456299098_3eaf656e12_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    low = low$V, g = g$V, m = m$V, gm = gm$V, gi = gi$V, 
    mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7135/7456300466_8393492196_o.png) 


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

![plot of chunk onerep](http://farm8.staticflickr.com/7258/7456400864_1748aa41da_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8018/7456453218_803de05dcd_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7276/7456482324_029edd79a8_o.png) 


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
<!-- Wed Jun 27 13:03:23 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 668.18 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 673.80 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 669.44 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 667.72 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 655.21 </TD> <TD align="right"> 668.76 </TD> <TD align="right"> 663.01 </TD> <TD align="right"> 666.28 </TD> <TD align="right"> 662.38 </TD> <TD align="right"> 668.38 </TD> <TD align="right"> 663.94 </TD> <TD align="right"> 665.90 </TD> <TD align="right"> 663.39 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 676.01 </TD> <TD align="right"> 669.55 </TD> <TD align="right"> 665.26 </TD> <TD align="right"> 669.83 </TD> <TD align="right"> 662.39 </TD> <TD align="right"> 670.77 </TD> <TD align="right"> 664.97 </TD> <TD align="right"> 671.02 </TD> <TD align="right"> 672.31 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 560.68 </TD> <TD align="right"> 561.81 </TD> <TD align="right"> 555.41 </TD> <TD align="right"> 586.41 </TD> <TD align="right"> 568.62 </TD> <TD align="right"> 585.38 </TD> <TD align="right"> 534.25 </TD> <TD align="right"> 583.42 </TD> <TD align="right"> 583.70 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 644.45 </TD> <TD align="right"> 651.91 </TD> <TD align="right"> 653.39 </TD> <TD align="right"> 650.76 </TD> <TD align="right"> 651.35 </TD> <TD align="right"> 651.67 </TD> <TD align="right"> 646.09 </TD> <TD align="right"> 647.85 </TD> <TD align="right"> 645.94 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 553.59 </TD> <TD align="right"> 553.69 </TD> <TD align="right"> 564.82 </TD> <TD align="right"> 550.48 </TD> <TD align="right"> 550.72 </TD> <TD align="right"> 566.30 </TD> <TD align="right"> 517.61 </TD> <TD align="right"> 587.03 </TD> <TD align="right"> 578.84 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 656.00 </TD> <TD align="right"> 621.05 </TD> <TD align="right"> 635.86 </TD> <TD align="right"> 637.60 </TD> <TD align="right"> 645.65 </TD> <TD align="right"> 650.77 </TD> <TD align="right"> 631.54 </TD> <TD align="right"> 651.00 </TD> <TD align="right"> 630.45 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 422.04 </TD> <TD align="right"> 402.17 </TD> <TD align="right"> 404.81 </TD> <TD align="right"> 454.65 </TD> <TD align="right"> 384.63 </TD> <TD align="right"> 457.73 </TD> <TD align="right"> 388.62 </TD> <TD align="right"> 461.04 </TD> <TD align="right"> 400.18 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 448.46 </TD> <TD align="right"> 390.15 </TD> <TD align="right"> 428.64 </TD> <TD align="right"> 460.10 </TD> <TD align="right"> 389.15 </TD> <TD align="right"> 452.20 </TD> <TD align="right"> 373.66 </TD> <TD align="right"> 442.54 </TD> <TD align="right"> 439.77 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), 
    dimnames = list(uncertainties, uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Wed Jun 27 13:03:23 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 23.05 </TD> <TD align="right"> 22.01 </TD> <TD align="right"> 20.66 </TD> <TD align="right"> 20.25 </TD> <TD align="right"> 19.36 </TD> <TD align="right"> 19.34 </TD> <TD align="right"> 20.80 </TD> <TD align="right"> 23.50 </TD> <TD align="right"> 20.45 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 103.22 </TD> <TD align="right"> 104.97 </TD> <TD align="right"> 100.88 </TD> <TD align="right"> 103.54 </TD> <TD align="right"> 97.84 </TD> <TD align="right"> 97.57 </TD> <TD align="right"> 95.14 </TD> <TD align="right"> 97.00 </TD> <TD align="right"> 95.69 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 53.07 </TD> <TD align="right"> 30.75 </TD> <TD align="right"> 85.80 </TD> <TD align="right"> 22.57 </TD> <TD align="right"> 31.87 </TD> <TD align="right"> 25.19 </TD> <TD align="right"> 109.12 </TD> <TD align="right"> 21.92 </TD> <TD align="right"> 19.33 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 13.13 </TD> <TD align="right"> 12.47 </TD> <TD align="right"> 11.41 </TD> <TD align="right"> 12.11 </TD> <TD align="right"> 12.56 </TD> <TD align="right"> 13.79 </TD> <TD align="right"> 12.82 </TD> <TD align="right"> 13.73 </TD> <TD align="right"> 12.49 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 103.31 </TD> <TD align="right"> 94.24 </TD> <TD align="right"> 105.53 </TD> <TD align="right"> 95.38 </TD> <TD align="right"> 88.98 </TD> <TD align="right"> 97.09 </TD> <TD align="right"> 138.85 </TD> <TD align="right"> 80.99 </TD> <TD align="right"> 95.29 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 93.57 </TD> <TD align="right"> 90.14 </TD> <TD align="right"> 98.76 </TD> <TD align="right"> 96.50 </TD> <TD align="right"> 110.70 </TD> <TD align="right"> 108.83 </TD> <TD align="right"> 87.15 </TD> <TD align="right"> 88.78 </TD> <TD align="right"> 93.71 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 170.63 </TD> <TD align="right"> 186.75 </TD> <TD align="right"> 170.15 </TD> <TD align="right"> 160.77 </TD> <TD align="right"> 184.13 </TD> <TD align="right"> 176.98 </TD> <TD align="right"> 188.16 </TD> <TD align="right"> 155.68 </TD> <TD align="right"> 186.02 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 149.85 </TD> <TD align="right"> 184.20 </TD> <TD align="right"> 175.38 </TD> <TD align="right"> 169.01 </TD> <TD align="right"> 172.25 </TD> <TD align="right"> 165.63 </TD> <TD align="right"> 178.70 </TD> <TD align="right"> 174.42 </TD> <TD align="right"> 180.79 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


