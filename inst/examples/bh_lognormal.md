






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

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, stock - x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7116/7463588954_baf23280c6_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, x_grid[value], color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8013/7463589282_56899d96b0_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7262/7463589528_d221a029bb_o.png) 


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

![plot of chunk onerep](http://farm9.staticflickr.com/8153/7463596070_0811da4dcb_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8167/7463598666_0713f06f6b_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8156/7463600340_70b4b5b51d_o.png) 


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
<!-- Thu Jun 28 16:24:10 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 318.18 </TD> <TD align="right"> 331.49 </TD> <TD align="right"> 331.73 </TD> <TD align="right"> 309.09 </TD> <TD align="right"> 331.80 </TD> <TD align="right"> 317.63 </TD> <TD align="right"> 324.24 </TD> <TD align="right"> 318.96 </TD> <TD align="right"> 319.21 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 315.26 </TD> <TD align="right"> 328.47 </TD> <TD align="right"> 328.25 </TD> <TD align="right"> 317.61 </TD> <TD align="right"> 333.27 </TD> <TD align="right"> 327.02 </TD> <TD align="right"> 318.65 </TD> <TD align="right"> 326.91 </TD> <TD align="right"> 323.43 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 489.02 </TD> <TD align="right"> 451.76 </TD> <TD align="right"> 461.19 </TD> <TD align="right"> 406.91 </TD> <TD align="right"> 440.00 </TD> <TD align="right"> 393.64 </TD> <TD align="right"> 444.60 </TD> <TD align="right"> 421.21 </TD> <TD align="right"> 416.46 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 200.13 </TD> <TD align="right"> 195.60 </TD> <TD align="right"> 159.62 </TD> <TD align="right"> 264.33 </TD> <TD align="right"> 148.01 </TD> <TD align="right"> 245.58 </TD> <TD align="right"> 158.74 </TD> <TD align="right"> 254.71 </TD> <TD align="right"> 233.83 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 288.01 </TD> <TD align="right"> 291.49 </TD> <TD align="right"> 280.08 </TD> <TD align="right"> 308.85 </TD> <TD align="right"> 283.60 </TD> <TD align="right"> 305.13 </TD> <TD align="right"> 269.88 </TD> <TD align="right"> 309.42 </TD> <TD align="right"> 305.92 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 316.07 </TD> <TD align="right"> 241.36 </TD> <TD align="right"> 201.05 </TD> <TD align="right"> 332.93 </TD> <TD align="right"> 199.35 </TD> <TD align="right"> 345.12 </TD> <TD align="right"> 167.60 </TD> <TD align="right"> 320.32 </TD> <TD align="right"> 325.72 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 385.44 </TD> <TD align="right"> 302.00 </TD> <TD align="right"> 310.08 </TD> <TD align="right"> 405.32 </TD> <TD align="right"> 301.00 </TD> <TD align="right"> 396.23 </TD> <TD align="right"> 283.87 </TD> <TD align="right"> 405.13 </TD> <TD align="right"> 401.91 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 176.47 </TD> <TD align="right"> 155.79 </TD> <TD align="right"> 165.19 </TD> <TD align="right"> 217.49 </TD> <TD align="right"> 155.83 </TD> <TD align="right"> 215.70 </TD> <TD align="right"> 151.01 </TD> <TD align="right"> 207.82 </TD> <TD align="right"> 204.85 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 238.84 </TD> <TD align="right"> 230.40 </TD> <TD align="right"> 185.19 </TD> <TD align="right"> 292.89 </TD> <TD align="right"> 216.60 </TD> <TD align="right"> 294.51 </TD> <TD align="right"> 189.44 </TD> <TD align="right"> 298.38 </TD> <TD align="right"> 236.78 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 28 16:24:10 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 33.67 </TD> <TD align="right"> 28.55 </TD> <TD align="right"> 37.90 </TD> <TD align="right"> 26.90 </TD> <TD align="right"> 26.20 </TD> <TD align="right"> 29.92 </TD> <TD align="right"> 22.21 </TD> <TD align="right"> 26.49 </TD> <TD align="right"> 26.57 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 169.73 </TD> <TD align="right"> 156.71 </TD> <TD align="right"> 159.37 </TD> <TD align="right"> 141.21 </TD> <TD align="right"> 151.03 </TD> <TD align="right"> 125.92 </TD> <TD align="right"> 125.84 </TD> <TD align="right"> 146.94 </TD> <TD align="right"> 139.48 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 60.74 </TD> <TD align="right"> 61.09 </TD> <TD align="right"> 67.94 </TD> <TD align="right"> 37.42 </TD> <TD align="right"> 55.36 </TD> <TD align="right"> 51.72 </TD> <TD align="right"> 63.99 </TD> <TD align="right"> 47.38 </TD> <TD align="right"> 58.45 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 46.98 </TD> <TD align="right"> 69.07 </TD> <TD align="right"> 81.77 </TD> <TD align="right"> 22.05 </TD> <TD align="right"> 73.43 </TD> <TD align="right"> 41.43 </TD> <TD align="right"> 76.54 </TD> <TD align="right"> 31.61 </TD> <TD align="right"> 36.41 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 158.87 </TD> <TD align="right"> 128.05 </TD> <TD align="right"> 127.50 </TD> <TD align="right"> 136.65 </TD> <TD align="right"> 138.57 </TD> <TD align="right"> 143.10 </TD> <TD align="right"> 103.42 </TD> <TD align="right"> 129.57 </TD> <TD align="right"> 124.54 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 162.28 </TD> <TD align="right"> 131.16 </TD> <TD align="right"> 150.05 </TD> <TD align="right"> 159.96 </TD> <TD align="right"> 155.50 </TD> <TD align="right"> 137.07 </TD> <TD align="right"> 121.44 </TD> <TD align="right"> 145.30 </TD> <TD align="right"> 155.75 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 64.12 </TD> <TD align="right"> 58.05 </TD> <TD align="right"> 67.45 </TD> <TD align="right"> 68.60 </TD> <TD align="right"> 60.55 </TD> <TD align="right"> 66.46 </TD> <TD align="right"> 54.65 </TD> <TD align="right"> 68.52 </TD> <TD align="right"> 64.94 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 131.60 </TD> <TD align="right"> 141.12 </TD> <TD align="right"> 109.00 </TD> <TD align="right"> 143.23 </TD> <TD align="right"> 121.42 </TD> <TD align="right"> 153.92 </TD> <TD align="right"> 111.48 </TD> <TD align="right"> 138.94 </TD> <TD align="right"> 113.67 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


