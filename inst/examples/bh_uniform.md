






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

![plot of chunk sethiplots](http://farm9.staticflickr.com/8002/7478125636_cb710ee9b4_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable), 
    shape = "+") + stat_smooth(aes(stock, x_grid[value], color = variable), 
    degree = 1, se = FALSE, span = 0.3) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7251/7478126126_3ba019a43c_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable), shape = "+") + 
    # stat_smooth(aes(stock, value, color=variable), degree=0, se=FALSE,
# span=0.15) +
ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7262/7478126470_d60c2f46d9_o.png) 








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

![plot of chunk onerep](http://farm9.staticflickr.com/8155/7478133886_57f03f4f27_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8142/7478137216_395a194a90_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8004/7478138736_5b4c1aeb1c_o.png) 


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
<!-- Sun Jul  1 03:55:12 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 31.61 </TD> <TD align="right"> 33.12 </TD> <TD align="right"> 32.96 </TD> <TD align="right"> 32.03 </TD> <TD align="right"> 32.98 </TD> <TD align="right"> 31.99 </TD> <TD align="right"> 32.68 </TD> <TD align="right"> 32.00 </TD> <TD align="right"> 31.76 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 31.68 </TD> <TD align="right"> 32.56 </TD> <TD align="right"> 32.90 </TD> <TD align="right"> 32.31 </TD> <TD align="right"> 32.89 </TD> <TD align="right"> 32.11 </TD> <TD align="right"> 32.26 </TD> <TD align="right"> 31.87 </TD> <TD align="right"> 32.18 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 31.70 </TD> <TD align="right"> 33.20 </TD> <TD align="right"> 32.67 </TD> <TD align="right"> 32.21 </TD> <TD align="right"> 32.91 </TD> <TD align="right"> 32.27 </TD> <TD align="right"> 32.75 </TD> <TD align="right"> 32.45 </TD> <TD align="right"> 30.41 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 22.24 </TD> <TD align="right"> 29.52 </TD> <TD align="right"> 30.06 </TD> <TD align="right"> 30.79 </TD> <TD align="right"> 29.04 </TD> <TD align="right"> 30.78 </TD> <TD align="right"> 30.85 </TD> <TD align="right"> 30.83 </TD> <TD align="right"> 30.71 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 28.46 </TD> <TD align="right"> 32.15 </TD> <TD align="right"> 32.06 </TD> <TD align="right"> 31.99 </TD> <TD align="right"> 32.14 </TD> <TD align="right"> 31.98 </TD> <TD align="right"> 31.88 </TD> <TD align="right"> 31.92 </TD> <TD align="right"> 31.85 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 20.96 </TD> <TD align="right"> 26.22 </TD> <TD align="right"> 28.06 </TD> <TD align="right"> 29.44 </TD> <TD align="right"> 28.97 </TD> <TD align="right"> 30.01 </TD> <TD align="right"> 29.79 </TD> <TD align="right"> 29.61 </TD> <TD align="right"> 29.32 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 27.72 </TD> <TD align="right"> 31.16 </TD> <TD align="right"> 30.84 </TD> <TD align="right"> 31.40 </TD> <TD align="right"> 31.03 </TD> <TD align="right"> 31.62 </TD> <TD align="right"> 31.89 </TD> <TD align="right"> 31.44 </TD> <TD align="right"> 31.44 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 18.86 </TD> <TD align="right"> 24.64 </TD> <TD align="right"> 27.04 </TD> <TD align="right"> 27.96 </TD> <TD align="right"> 26.00 </TD> <TD align="right"> 28.52 </TD> <TD align="right"> 26.73 </TD> <TD align="right"> 28.88 </TD> <TD align="right"> 27.52 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 19.22 </TD> <TD align="right"> 22.36 </TD> <TD align="right"> 24.71 </TD> <TD align="right"> 25.80 </TD> <TD align="right"> 24.77 </TD> <TD align="right"> 28.07 </TD> <TD align="right"> 26.58 </TD> <TD align="right"> 26.39 </TD> <TD align="right"> 27.04 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Sun Jul  1 03:55:12 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> gro </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 1.61 </TD> <TD align="right"> 1.56 </TD> <TD align="right"> 1.67 </TD> <TD align="right"> 1.47 </TD> <TD align="right"> 1.67 </TD> <TD align="right"> 1.61 </TD> <TD align="right"> 1.60 </TD> <TD align="right"> 1.46 </TD> <TD align="right"> 1.79 </TD> </TR>
  <TR> <TD align="right"> gro </TD> <TD align="right"> 7.36 </TD> <TD align="right"> 8.36 </TD> <TD align="right"> 8.36 </TD> <TD align="right"> 6.83 </TD> <TD align="right"> 6.71 </TD> <TD align="right"> 6.78 </TD> <TD align="right"> 7.17 </TD> <TD align="right"> 7.47 </TD> <TD align="right"> 7.56 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 6.96 </TD> <TD align="right"> 4.26 </TD> <TD align="right"> 3.23 </TD> <TD align="right"> 0.79 </TD> <TD align="right"> 5.35 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 1.08 </TD> <TD align="right"> 0.86 </TD> <TD align="right"> 0.90 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 5.48 </TD> <TD align="right"> 0.87 </TD> <TD align="right"> 0.93 </TD> <TD align="right"> 0.37 </TD> <TD align="right"> 0.89 </TD> <TD align="right"> 0.42 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 0.38 </TD> <TD align="right"> 0.33 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 7.84 </TD> <TD align="right"> 7.32 </TD> <TD align="right"> 6.86 </TD> <TD align="right"> 6.52 </TD> <TD align="right"> 8.29 </TD> <TD align="right"> 7.99 </TD> <TD align="right"> 6.89 </TD> <TD align="right"> 6.50 </TD> <TD align="right"> 6.24 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 8.13 </TD> <TD align="right"> 6.72 </TD> <TD align="right"> 7.86 </TD> <TD align="right"> 7.54 </TD> <TD align="right"> 7.92 </TD> <TD align="right"> 7.22 </TD> <TD align="right"> 8.49 </TD> <TD align="right"> 7.93 </TD> <TD align="right"> 8.62 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 7.27 </TD> <TD align="right"> 8.04 </TD> <TD align="right"> 6.31 </TD> <TD align="right"> 5.39 </TD> <TD align="right"> 7.68 </TD> <TD align="right"> 4.98 </TD> <TD align="right"> 7.12 </TD> <TD align="right"> 5.13 </TD> <TD align="right"> 6.20 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 8.02 </TD> <TD align="right"> 8.45 </TD> <TD align="right"> 9.64 </TD> <TD align="right"> 8.31 </TD> <TD align="right"> 8.95 </TD> <TD align="right"> 8.06 </TD> <TD align="right"> 8.96 </TD> <TD align="right"> 8.20 </TD> <TD align="right"> 7.83 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


