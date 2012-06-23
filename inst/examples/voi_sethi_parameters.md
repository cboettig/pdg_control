






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
det <- scenario(0.00101, 0.001, 0.001)
```

```
Library ggplot2 loaded.
```






```r
g <- scenario(lvl, 0.001, 0.001)
```

```
Library ggplot2 loaded.
```






```r
m <- scenario(0.001, lvl, 0.001)
```

```
Library ggplot2 loaded.
```






```r
i <- scenario(0.001, 0.001, lvl)
```

```
Library ggplot2 loaded.
```






```r
gm <- scenario(lvl, lvl, 0.001)
```

```
Library ggplot2 loaded.
```






```r
gi <- scenario(lvl, 0.001, lvl)
```

```
Library ggplot2 loaded.
```






```r
mi <- scenario(0.001, lvl, lvl)
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





### plots



```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, det = det$D[, 
    1], g = g$D[, 1], m = m$D[, 1], i = m$D[, 1], gm = gm$D[, 
    1], gi = gi$D[, 1], gmi = gmi$D[, 1]), id = "stock")

ggplot(policy) + geom_point(aes(stock, stock - 
    x_grid[value], color = variable)) + geom_smooth(aes(stock, 
    stock - x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7256/7424377052_9b7e722a9f_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], 
    color = variable)) + geom_smooth(aes(stock, x_grid[value], 
    color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7106/7424377306_667cca255c_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, gmi = gmi$V), 
    id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm6.staticflickr.com/5465/7424377532_31eaebd7ef_o.png) 


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
policyfn <- list(det = det, g = g, m = m, i = i, 
    gm = gm, gi = gi, mi = mi, gmi = gmi)
noise <- list(s0.001 = c(0.001, 0.001, 0.001), 
    sg = c(lvl, 0.001, 0.001), sm = c(0.001, lvl, 0.001), 
    si = c(0.001, 0.001, lvl), sgm = c(lvl, lvl, 0.001), 
    sgi = c(lvl, 0.001, lvl), smi = c(0.001, lvl, lvl), sgmi = c(lvl, 
        lvl, lvl))
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

![plot of chunk onerep](http://farm6.staticflickr.com/5470/7424381870_5cf372b85c_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7273/7424384394_623922d85c_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8151/7424385896_34e667446e_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by = uncertainty]
sds <- profits[, sd(V1), by = uncertainty]
```






```r
require(xtable)
uncertainties <- c("det", "growth", "measure", 
    "impl", "growth+measure", "growth+impl", "impl+measure", 
    "all")
print(xtable(matrix(means$V1, nrow = 8, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Sat Jun 23 12:55:12 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> growth </TH> <TH> measure </TH> <TH> impl </TH> <TH> growth+measure </TH> <TH> growth+impl </TH> <TH> impl+measure </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 668.16 </TD> <TD align="right"> 674.53 </TD> <TD align="right"> 674.17 </TD> <TD align="right"> 672.73 </TD> <TD align="right"> 673.59 </TD> <TD align="right"> 668.30 </TD> <TD align="right"> 673.29 </TD> <TD align="right"> 663.49 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 663.73 </TD> <TD align="right"> 691.08 </TD> <TD align="right"> 664.49 </TD> <TD align="right"> 684.95 </TD> <TD align="right"> 692.24 </TD> <TD align="right"> 682.32 </TD> <TD align="right"> 643.37 </TD> <TD align="right"> 658.22 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 556.42 </TD> <TD align="right"> 563.60 </TD> <TD align="right"> 586.39 </TD> <TD align="right"> 550.40 </TD> <TD align="right"> 592.05 </TD> <TD align="right"> 515.50 </TD> <TD align="right"> 580.60 </TD> <TD align="right"> 585.14 </TD> </TR>
  <TR> <TD align="right"> impl </TD> <TD align="right"> 644.11 </TD> <TD align="right"> 650.11 </TD> <TD align="right"> 651.09 </TD> <TD align="right"> 650.42 </TD> <TD align="right"> 649.24 </TD> <TD align="right"> 648.52 </TD> <TD align="right"> 648.75 </TD> <TD align="right"> 641.25 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 566.73 </TD> <TD align="right"> 575.88 </TD> <TD align="right"> 582.07 </TD> <TD align="right"> 557.54 </TD> <TD align="right"> 588.06 </TD> <TD align="right"> 539.90 </TD> <TD align="right"> 581.51 </TD> <TD align="right"> 564.03 </TD> </TR>
  <TR> <TD align="right"> growth+impl </TD> <TD align="right"> 634.21 </TD> <TD align="right"> 641.05 </TD> <TD align="right"> 624.10 </TD> <TD align="right"> 639.99 </TD> <TD align="right"> 640.94 </TD> <TD align="right"> 635.67 </TD> <TD align="right"> 647.07 </TD> <TD align="right"> 620.94 </TD> </TR>
  <TR> <TD align="right"> impl+measure </TD> <TD align="right"> 405.89 </TD> <TD align="right"> 389.87 </TD> <TD align="right"> 437.02 </TD> <TD align="right"> 389.76 </TD> <TD align="right"> 442.34 </TD> <TD align="right"> 419.33 </TD> <TD align="right"> 452.07 </TD> <TD align="right"> 438.89 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 423.59 </TD> <TD align="right"> 432.70 </TD> <TD align="right"> 470.25 </TD> <TD align="right"> 399.05 </TD> <TD align="right"> 449.59 </TD> <TD align="right"> 355.47 </TD> <TD align="right"> 458.69 </TD> <TD align="right"> 414.71 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = 8, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Sat Jun 23 12:55:12 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> growth </TH> <TH> measure </TH> <TH> impl </TH> <TH> growth+measure </TH> <TH> growth+impl </TH> <TH> impl+measure </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.08 </TD> <TD align="right"> 0.41 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.19 </TD> <TD align="right"> 0.46 </TD> <TD align="right"> 0.31 </TD> <TD align="right"> 0.47 </TD> <TD align="right"> 0.70 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 110.61 </TD> <TD align="right"> 102.34 </TD> <TD align="right"> 101.63 </TD> <TD align="right"> 101.61 </TD> <TD align="right"> 112.49 </TD> <TD align="right"> 106.85 </TD> <TD align="right"> 103.35 </TD> <TD align="right"> 82.97 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 34.03 </TD> <TD align="right"> 31.69 </TD> <TD align="right"> 20.59 </TD> <TD align="right"> 69.68 </TD> <TD align="right"> 23.47 </TD> <TD align="right"> 131.48 </TD> <TD align="right"> 24.15 </TD> <TD align="right"> 18.55 </TD> </TR>
  <TR> <TD align="right"> impl </TD> <TD align="right"> 13.26 </TD> <TD align="right"> 11.88 </TD> <TD align="right"> 13.04 </TD> <TD align="right"> 12.59 </TD> <TD align="right"> 12.99 </TD> <TD align="right"> 10.50 </TD> <TD align="right"> 13.87 </TD> <TD align="right"> 13.71 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 95.77 </TD> <TD align="right"> 93.08 </TD> <TD align="right"> 94.80 </TD> <TD align="right"> 93.17 </TD> <TD align="right"> 100.72 </TD> <TD align="right"> 114.70 </TD> <TD align="right"> 85.79 </TD> <TD align="right"> 81.50 </TD> </TR>
  <TR> <TD align="right"> growth+impl </TD> <TD align="right"> 107.34 </TD> <TD align="right"> 90.39 </TD> <TD align="right"> 100.70 </TD> <TD align="right"> 98.66 </TD> <TD align="right"> 97.27 </TD> <TD align="right"> 99.37 </TD> <TD align="right"> 110.67 </TD> <TD align="right"> 90.89 </TD> </TR>
  <TR> <TD align="right"> impl+measure </TD> <TD align="right"> 174.01 </TD> <TD align="right"> 185.83 </TD> <TD align="right"> 169.40 </TD> <TD align="right"> 181.17 </TD> <TD align="right"> 170.39 </TD> <TD align="right"> 169.55 </TD> <TD align="right"> 157.74 </TD> <TD align="right"> 179.14 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 166.27 </TD> <TD align="right"> 168.65 </TD> <TD align="right"> 143.18 </TD> <TD align="right"> 182.11 </TD> <TD align="right"> 157.54 </TD> <TD align="right"> 183.03 </TD> <TD align="right"> 173.67 </TD> <TD align="right"> 173.95 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


