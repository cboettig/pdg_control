






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
    
    z_g <- function() rlnorm(1, 0, policy_g)
    z_m <- function() rlnorm(1, 0, policy_m)
    z_i <- function() rlnorm(1, 0, policy_i)
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps = 20000)
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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7125/7469644504_59df8f5dce_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7140/7469644784_1159eb4c28_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7136/7469645060_01c550d491_o.png) 




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

![plot of chunk onerep](http://farm8.staticflickr.com/7132/7469650652_eccdb2d244_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7247/7469653076_d8a3b76e43_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7111/7469654324_5eff07b886_o.png) 


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
<!-- Fri Jun 29 17:13:24 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 318.18 </TD> <TD align="right"> 331.82 </TD> <TD align="right"> 329.84 </TD> <TD align="right"> 331.82 </TD> <TD align="right"> 331.56 </TD> <TD align="right"> 327.27 </TD> <TD align="right"> 326.85 </TD> <TD align="right"> 330.30 </TD> <TD align="right"> 328.64 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 315.26 </TD> <TD align="right"> 328.61 </TD> <TD align="right"> 328.10 </TD> <TD align="right"> 331.86 </TD> <TD align="right"> 334.49 </TD> <TD align="right"> 331.74 </TD> <TD align="right"> 323.59 </TD> <TD align="right"> 334.02 </TD> <TD align="right"> 328.43 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 351.77 </TD> <TD align="right"> 352.50 </TD> <TD align="right"> 357.22 </TD> <TD align="right"> 348.45 </TD> <TD align="right"> 350.74 </TD> <TD align="right"> 346.36 </TD> <TD align="right"> 363.11 </TD> <TD align="right"> 354.33 </TD> <TD align="right"> 350.50 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 304.16 </TD> <TD align="right"> 310.53 </TD> <TD align="right"> 308.81 </TD> <TD align="right"> 318.18 </TD> <TD align="right"> 308.91 </TD> <TD align="right"> 313.09 </TD> <TD align="right"> 307.81 </TD> <TD align="right"> 313.05 </TD> <TD align="right"> 314.81 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 309.39 </TD> <TD align="right"> 326.71 </TD> <TD align="right"> 326.07 </TD> <TD align="right"> 327.25 </TD> <TD align="right"> 326.60 </TD> <TD align="right"> 324.04 </TD> <TD align="right"> 322.71 </TD> <TD align="right"> 326.22 </TD> <TD align="right"> 324.88 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 337.97 </TD> <TD align="right"> 337.80 </TD> <TD align="right"> 328.30 </TD> <TD align="right"> 339.78 </TD> <TD align="right"> 307.55 </TD> <TD align="right"> 333.02 </TD> <TD align="right"> 324.93 </TD> <TD align="right"> 331.86 </TD> <TD align="right"> 342.77 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 343.60 </TD> <TD align="right"> 343.46 </TD> <TD align="right"> 349.40 </TD> <TD align="right"> 354.16 </TD> <TD align="right"> 346.61 </TD> <TD align="right"> 347.43 </TD> <TD align="right"> 345.61 </TD> <TD align="right"> 353.86 </TD> <TD align="right"> 360.82 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 293.85 </TD> <TD align="right"> 284.07 </TD> <TD align="right"> 294.72 </TD> <TD align="right"> 295.79 </TD> <TD align="right"> 290.21 </TD> <TD align="right"> 307.00 </TD> <TD align="right"> 296.86 </TD> <TD align="right"> 311.70 </TD> <TD align="right"> 308.94 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 325.79 </TD> <TD align="right"> 318.84 </TD> <TD align="right"> 303.77 </TD> <TD align="right"> 329.83 </TD> <TD align="right"> 323.60 </TD> <TD align="right"> 329.15 </TD> <TD align="right"> 316.73 </TD> <TD align="right"> 335.38 </TD> <TD align="right"> 317.80 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Fri Jun 29 17:13:24 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 33.67 </TD> <TD align="right"> 28.31 </TD> <TD align="right"> 39.45 </TD> <TD align="right"> 25.76 </TD> <TD align="right"> 27.40 </TD> <TD align="right"> 31.41 </TD> <TD align="right"> 27.39 </TD> <TD align="right"> 27.38 </TD> <TD align="right"> 30.50 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 63.60 </TD> <TD align="right"> 58.09 </TD> <TD align="right"> 60.52 </TD> <TD align="right"> 51.68 </TD> <TD align="right"> 54.95 </TD> <TD align="right"> 57.92 </TD> <TD align="right"> 56.75 </TD> <TD align="right"> 56.44 </TD> <TD align="right"> 57.15 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 22.49 </TD> <TD align="right"> 31.81 </TD> <TD align="right"> 37.86 </TD> <TD align="right"> 13.33 </TD> <TD align="right"> 29.57 </TD> <TD align="right"> 22.47 </TD> <TD align="right"> 36.29 </TD> <TD align="right"> 23.35 </TD> <TD align="right"> 9.06 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 8.00 </TD> <TD align="right"> 6.10 </TD> <TD align="right"> 4.74 </TD> <TD align="right"> 5.88 </TD> <TD align="right"> 6.08 </TD> <TD align="right"> 5.93 </TD> <TD align="right"> 6.36 </TD> <TD align="right"> 6.67 </TD> <TD align="right"> 5.67 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 61.95 </TD> <TD align="right"> 60.45 </TD> <TD align="right"> 71.49 </TD> <TD align="right"> 52.83 </TD> <TD align="right"> 83.75 </TD> <TD align="right"> 59.95 </TD> <TD align="right"> 70.35 </TD> <TD align="right"> 51.09 </TD> <TD align="right"> 57.76 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 57.71 </TD> <TD align="right"> 55.48 </TD> <TD align="right"> 52.68 </TD> <TD align="right"> 53.45 </TD> <TD align="right"> 58.05 </TD> <TD align="right"> 56.69 </TD> <TD align="right"> 56.87 </TD> <TD align="right"> 53.75 </TD> <TD align="right"> 61.30 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 41.20 </TD> <TD align="right"> 65.59 </TD> <TD align="right"> 59.68 </TD> <TD align="right"> 52.05 </TD> <TD align="right"> 57.75 </TD> <TD align="right"> 29.71 </TD> <TD align="right"> 49.17 </TD> <TD align="right"> 15.12 </TD> <TD align="right"> 23.80 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 62.71 </TD> <TD align="right"> 77.63 </TD> <TD align="right"> 84.30 </TD> <TD align="right"> 61.97 </TD> <TD align="right"> 60.77 </TD> <TD align="right"> 58.59 </TD> <TD align="right"> 69.11 </TD> <TD align="right"> 56.58 </TD> <TD align="right"> 69.44 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


