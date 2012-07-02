






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

![plot of chunk sethiplots](http://farm9.staticflickr.com/8007/7475801558_15bc551873_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7250/7475802062_61af36a4d7_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8020/7475802392_b76d7bab2a_o.png) 




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

![plot of chunk onerep](http://farm8.staticflickr.com/7280/7475809270_d5055d7799_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8166/7475812330_93fe57723e_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8006/7475813770_314e6c1976_o.png) 


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
<!-- Sat Jun 30 17:57:42 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 315.11 </TD> <TD align="right"> 331.21 </TD> <TD align="right"> 330.25 </TD> <TD align="right"> 331.76 </TD> <TD align="right"> 329.96 </TD> <TD align="right"> 328.74 </TD> <TD align="right"> 325.71 </TD> <TD align="right"> 329.74 </TD> <TD align="right"> 329.93 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 311.99 </TD> <TD align="right"> 328.73 </TD> <TD align="right"> 332.36 </TD> <TD align="right"> 330.70 </TD> <TD align="right"> 328.11 </TD> <TD align="right"> 324.75 </TD> <TD align="right"> 318.74 </TD> <TD align="right"> 327.90 </TD> <TD align="right"> 326.07 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 344.73 </TD> <TD align="right"> 357.72 </TD> <TD align="right"> 352.08 </TD> <TD align="right"> 350.68 </TD> <TD align="right"> 350.97 </TD> <TD align="right"> 359.95 </TD> <TD align="right"> 365.24 </TD> <TD align="right"> 351.98 </TD> <TD align="right"> 360.70 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 240.74 </TD> <TD align="right"> 314.51 </TD> <TD align="right"> 311.34 </TD> <TD align="right"> 314.71 </TD> <TD align="right"> 302.35 </TD> <TD align="right"> 316.16 </TD> <TD align="right"> 312.73 </TD> <TD align="right"> 314.68 </TD> <TD align="right"> 315.62 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 303.84 </TD> <TD align="right"> 327.59 </TD> <TD align="right"> 325.44 </TD> <TD align="right"> 326.85 </TD> <TD align="right"> 325.47 </TD> <TD align="right"> 324.77 </TD> <TD align="right"> 321.94 </TD> <TD align="right"> 325.87 </TD> <TD align="right"> 323.26 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 270.98 </TD> <TD align="right"> 320.56 </TD> <TD align="right"> 333.18 </TD> <TD align="right"> 334.54 </TD> <TD align="right"> 329.80 </TD> <TD align="right"> 329.18 </TD> <TD align="right"> 330.74 </TD> <TD align="right"> 333.55 </TD> <TD align="right"> 332.63 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 327.57 </TD> <TD align="right"> 351.12 </TD> <TD align="right"> 355.88 </TD> <TD align="right"> 337.88 </TD> <TD align="right"> 343.83 </TD> <TD align="right"> 345.44 </TD> <TD align="right"> 344.60 </TD> <TD align="right"> 344.93 </TD> <TD align="right"> 346.90 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 219.49 </TD> <TD align="right"> 299.67 </TD> <TD align="right"> 296.66 </TD> <TD align="right"> 306.16 </TD> <TD align="right"> 290.73 </TD> <TD align="right"> 310.83 </TD> <TD align="right"> 294.66 </TD> <TD align="right"> 307.54 </TD> <TD align="right"> 302.51 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 247.27 </TD> <TD align="right"> 302.99 </TD> <TD align="right"> 305.50 </TD> <TD align="right"> 340.55 </TD> <TD align="right"> 316.67 </TD> <TD align="right"> 329.26 </TD> <TD align="right"> 316.54 </TD> <TD align="right"> 318.11 </TD> <TD align="right"> 318.05 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Sat Jun 30 17:57:42 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 31.46 </TD> <TD align="right"> 26.27 </TD> <TD align="right"> 25.60 </TD> <TD align="right"> 25.52 </TD> <TD align="right"> 27.63 </TD> <TD align="right"> 29.48 </TD> <TD align="right"> 28.63 </TD> <TD align="right"> 27.36 </TD> <TD align="right"> 32.19 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 54.40 </TD> <TD align="right"> 58.85 </TD> <TD align="right"> 58.67 </TD> <TD align="right"> 56.76 </TD> <TD align="right"> 54.23 </TD> <TD align="right"> 57.81 </TD> <TD align="right"> 65.36 </TD> <TD align="right"> 63.00 </TD> <TD align="right"> 56.55 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 75.49 </TD> <TD align="right"> 10.96 </TD> <TD align="right"> 26.40 </TD> <TD align="right"> 14.21 </TD> <TD align="right"> 50.42 </TD> <TD align="right"> 12.70 </TD> <TD align="right"> 13.23 </TD> <TD align="right"> 18.50 </TD> <TD align="right"> 7.92 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 37.14 </TD> <TD align="right"> 5.10 </TD> <TD align="right"> 6.24 </TD> <TD align="right"> 5.83 </TD> <TD align="right"> 6.04 </TD> <TD align="right"> 6.10 </TD> <TD align="right"> 6.58 </TD> <TD align="right"> 5.70 </TD> <TD align="right"> 6.43 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 88.37 </TD> <TD align="right"> 57.82 </TD> <TD align="right"> 74.88 </TD> <TD align="right"> 56.70 </TD> <TD align="right"> 69.83 </TD> <TD align="right"> 58.81 </TD> <TD align="right"> 58.99 </TD> <TD align="right"> 51.95 </TD> <TD align="right"> 59.98 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 67.55 </TD> <TD align="right"> 53.10 </TD> <TD align="right"> 57.78 </TD> <TD align="right"> 56.65 </TD> <TD align="right"> 50.32 </TD> <TD align="right"> 59.10 </TD> <TD align="right"> 51.50 </TD> <TD align="right"> 56.36 </TD> <TD align="right"> 52.77 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 77.98 </TD> <TD align="right"> 45.23 </TD> <TD align="right"> 48.01 </TD> <TD align="right"> 32.03 </TD> <TD align="right"> 57.81 </TD> <TD align="right"> 30.14 </TD> <TD align="right"> 50.16 </TD> <TD align="right"> 28.78 </TD> <TD align="right"> 37.60 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 79.56 </TD> <TD align="right"> 73.74 </TD> <TD align="right"> 76.26 </TD> <TD align="right"> 58.01 </TD> <TD align="right"> 67.29 </TD> <TD align="right"> 65.02 </TD> <TD align="right"> 61.75 </TD> <TD align="right"> 64.58 </TD> <TD align="right"> 62.05 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


