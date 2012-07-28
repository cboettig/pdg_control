






# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   




Chose the state equation / population dynamics function



```r
f <- BevHolt
```





With parameters `A` = `1.5` and `B` = `0.05`.



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




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `15` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7249/7498037722_87631d7cc6_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8282/7498038144_e671abf197_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8144/7498038410_c1594ace34_o.png) 








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

![plot of chunk onerep](http://farm8.staticflickr.com/7128/7498044450_38f7d32d7c_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8291/7498047150_f355f4c6ed_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8152/7498048494_f310d21c27_o.png) 


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
<!-- Tue Jul  3 17:42:04 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 31.24 </TD> <TD align="right"> 32.63 </TD> <TD align="right"> 32.36 </TD> <TD align="right"> 31.57 </TD> <TD align="right"> 32.59 </TD> <TD align="right"> 31.46 </TD> <TD align="right"> 31.94 </TD> <TD align="right"> 31.16 </TD> <TD align="right"> 31.19 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 31.11 </TD> <TD align="right"> 32.42 </TD> <TD align="right"> 31.81 </TD> <TD align="right"> 31.66 </TD> <TD align="right"> 32.27 </TD> <TD align="right"> 31.29 </TD> <TD align="right"> 31.74 </TD> <TD align="right"> 31.34 </TD> <TD align="right"> 31.61 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 29.99 </TD> <TD align="right"> 33.21 </TD> <TD align="right"> 31.82 </TD> <TD align="right"> 30.24 </TD> <TD align="right"> 32.06 </TD> <TD align="right"> 30.20 </TD> <TD align="right"> 31.76 </TD> <TD align="right"> 30.90 </TD> <TD align="right"> 31.30 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 20.60 </TD> <TD align="right"> 27.42 </TD> <TD align="right"> 29.58 </TD> <TD align="right"> 30.08 </TD> <TD align="right"> 29.76 </TD> <TD align="right"> 30.58 </TD> <TD align="right"> 30.19 </TD> <TD align="right"> 30.38 </TD> <TD align="right"> 30.20 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 29.17 </TD> <TD align="right"> 31.74 </TD> <TD align="right"> 31.55 </TD> <TD align="right"> 31.50 </TD> <TD align="right"> 31.83 </TD> <TD align="right"> 31.42 </TD> <TD align="right"> 31.37 </TD> <TD align="right"> 31.24 </TD> <TD align="right"> 31.21 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 20.61 </TD> <TD align="right"> 27.46 </TD> <TD align="right"> 28.47 </TD> <TD align="right"> 28.71 </TD> <TD align="right"> 28.60 </TD> <TD align="right"> 29.39 </TD> <TD align="right"> 29.94 </TD> <TD align="right"> 30.91 </TD> <TD align="right"> 28.83 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 26.13 </TD> <TD align="right"> 30.89 </TD> <TD align="right"> 30.95 </TD> <TD align="right"> 30.13 </TD> <TD align="right"> 30.13 </TD> <TD align="right"> 30.31 </TD> <TD align="right"> 31.08 </TD> <TD align="right"> 31.44 </TD> <TD align="right"> 30.14 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 18.57 </TD> <TD align="right"> 23.71 </TD> <TD align="right"> 25.29 </TD> <TD align="right"> 27.16 </TD> <TD align="right"> 26.28 </TD> <TD align="right"> 27.98 </TD> <TD align="right"> 25.70 </TD> <TD align="right"> 28.23 </TD> <TD align="right"> 27.30 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 17.94 </TD> <TD align="right"> 22.77 </TD> <TD align="right"> 24.74 </TD> <TD align="right"> 26.84 </TD> <TD align="right"> 25.91 </TD> <TD align="right"> 26.70 </TD> <TD align="right"> 24.17 </TD> <TD align="right"> 26.33 </TD> <TD align="right"> 28.13 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Tue Jul  3 17:42:04 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 1.55 </TD> <TD align="right"> 1.47 </TD> <TD align="right"> 1.59 </TD> <TD align="right"> 1.50 </TD> <TD align="right"> 1.69 </TD> <TD align="right"> 1.64 </TD> <TD align="right"> 1.96 </TD> <TD align="right"> 1.80 </TD> <TD align="right"> 1.85 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 7.49 </TD> <TD align="right"> 8.72 </TD> <TD align="right"> 6.54 </TD> <TD align="right"> 7.44 </TD> <TD align="right"> 6.67 </TD> <TD align="right"> 7.37 </TD> <TD align="right"> 8.42 </TD> <TD align="right"> 8.46 </TD> <TD align="right"> 7.53 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 7.59 </TD> <TD align="right"> 6.60 </TD> <TD align="right"> 3.06 </TD> <TD align="right"> 0.89 </TD> <TD align="right"> 2.38 </TD> <TD align="right"> 0.80 </TD> <TD align="right"> 1.29 </TD> <TD align="right"> 0.88 </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 4.20 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 0.98 </TD> <TD align="right"> 0.29 </TD> <TD align="right"> 0.79 </TD> <TD align="right"> 0.39 </TD> <TD align="right"> 0.96 </TD> <TD align="right"> 0.40 </TD> <TD align="right"> 0.53 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 8.09 </TD> <TD align="right"> 7.17 </TD> <TD align="right"> 6.51 </TD> <TD align="right"> 5.67 </TD> <TD align="right"> 5.77 </TD> <TD align="right"> 6.24 </TD> <TD align="right"> 7.32 </TD> <TD align="right"> 7.71 </TD> <TD align="right"> 6.64 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 7.69 </TD> <TD align="right"> 7.82 </TD> <TD align="right"> 7.28 </TD> <TD align="right"> 6.26 </TD> <TD align="right"> 7.43 </TD> <TD align="right"> 8.15 </TD> <TD align="right"> 8.27 </TD> <TD align="right"> 7.81 </TD> <TD align="right"> 8.31 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 7.08 </TD> <TD align="right"> 8.29 </TD> <TD align="right"> 7.25 </TD> <TD align="right"> 5.75 </TD> <TD align="right"> 6.28 </TD> <TD align="right"> 5.20 </TD> <TD align="right"> 7.31 </TD> <TD align="right"> 4.85 </TD> <TD align="right"> 5.97 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 8.09 </TD> <TD align="right"> 8.48 </TD> <TD align="right"> 8.50 </TD> <TD align="right"> 8.45 </TD> <TD align="right"> 8.22 </TD> <TD align="right"> 7.34 </TD> <TD align="right"> 8.94 </TD> <TD align="right"> 8.05 </TD> <TD align="right"> 8.01 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


