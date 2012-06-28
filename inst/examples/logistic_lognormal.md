






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

![plot of chunk sethiplots](http://farm9.staticflickr.com/8016/7460761844_d21ea3b767_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, x_grid[value], color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8026/7460762430_75f0c01bf6_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8008/7460762810_4107b372dc_o.png) 


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

![plot of chunk onerep](http://farm8.staticflickr.com/7137/7460771974_de788a00cb_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8006/7460775812_8d97cf8c57_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7247/7460777848_059dec5c73_o.png) 


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
<!-- Thu Jun 28 07:16:16 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 668.18 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 648.40 </TD> <TD align="right"> 602.80 </TD> <TD align="right"> 650.97 </TD> <TD align="right"> 620.99 </TD> <TD align="right"> 601.50 </TD> <TD align="right"> 629.59 </TD> <TD align="right"> 636.14 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 653.97 </TD> <TD align="right"> 661.49 </TD> <TD align="right"> 640.73 </TD> <TD align="right"> 619.19 </TD> <TD align="right"> 631.93 </TD> <TD align="right"> 635.22 </TD> <TD align="right"> 601.78 </TD> <TD align="right"> 633.96 </TD> <TD align="right"> 637.98 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 738.01 </TD> <TD align="right"> 774.61 </TD> <TD align="right"> 833.42 </TD> <TD align="right"> 529.27 </TD> <TD align="right"> 793.47 </TD> <TD align="right"> 551.40 </TD> <TD align="right"> 744.70 </TD> <TD align="right"> 532.50 </TD> <TD align="right"> 523.20 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 330.53 </TD> <TD align="right"> 261.09 </TD> <TD align="right"> 194.14 </TD> <TD align="right"> 514.62 </TD> <TD align="right"> 190.73 </TD> <TD align="right"> 522.37 </TD> <TD align="right"> 178.77 </TD> <TD align="right"> 518.50 </TD> <TD align="right"> 486.59 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 509.36 </TD> <TD align="right"> 482.41 </TD> <TD align="right"> 382.99 </TD> <TD align="right"> 571.14 </TD> <TD align="right"> 413.74 </TD> <TD align="right"> 583.63 </TD> <TD align="right"> 350.74 </TD> <TD align="right"> 586.21 </TD> <TD align="right"> 581.64 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 374.12 </TD> <TD align="right"> 357.62 </TD> <TD align="right"> 235.62 </TD> <TD align="right"> 424.01 </TD> <TD align="right"> 266.18 </TD> <TD align="right"> 425.98 </TD> <TD align="right"> 230.45 </TD> <TD align="right"> 445.87 </TD> <TD align="right"> 455.21 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 445.95 </TD> <TD align="right"> 446.30 </TD> <TD align="right"> 384.13 </TD> <TD align="right"> 539.77 </TD> <TD align="right"> 411.54 </TD> <TD align="right"> 539.21 </TD> <TD align="right"> 307.56 </TD> <TD align="right"> 498.22 </TD> <TD align="right"> 508.85 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 219.75 </TD> <TD align="right"> 195.19 </TD> <TD align="right"> 185.59 </TD> <TD align="right"> 394.42 </TD> <TD align="right"> 187.69 </TD> <TD align="right"> 370.92 </TD> <TD align="right"> 173.04 </TD> <TD align="right"> 369.77 </TD> <TD align="right"> 365.00 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 253.94 </TD> <TD align="right"> 245.79 </TD> <TD align="right"> 230.41 </TD> <TD align="right"> 347.08 </TD> <TD align="right"> 207.10 </TD> <TD align="right"> 342.81 </TD> <TD align="right"> 207.81 </TD> <TD align="right"> 347.40 </TD> <TD align="right"> 374.76 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 28 07:16:16 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth_measure </TH> <TH> growth_implement </TH> <TH> measure_implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 36.76 </TD> <TD align="right"> 36.30 </TD> <TD align="right"> 34.41 </TD> <TD align="right"> 34.54 </TD> <TD align="right"> 79.82 </TD> <TD align="right"> 37.10 </TD> <TD align="right"> 29.33 </TD> <TD align="right"> 37.19 </TD> <TD align="right"> 34.00 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 278.45 </TD> <TD align="right"> 197.65 </TD> <TD align="right"> 185.73 </TD> <TD align="right"> 282.77 </TD> <TD align="right"> 184.89 </TD> <TD align="right"> 292.57 </TD> <TD align="right"> 150.30 </TD> <TD align="right"> 315.74 </TD> <TD align="right"> 280.83 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 129.66 </TD> <TD align="right"> 101.53 </TD> <TD align="right"> 114.89 </TD> <TD align="right"> 61.84 </TD> <TD align="right"> 107.24 </TD> <TD align="right"> 75.87 </TD> <TD align="right"> 94.16 </TD> <TD align="right"> 73.02 </TD> <TD align="right"> 95.13 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 186.39 </TD> <TD align="right"> 190.43 </TD> <TD align="right"> 196.93 </TD> <TD align="right"> 111.58 </TD> <TD align="right"> 197.95 </TD> <TD align="right"> 73.62 </TD> <TD align="right"> 193.47 </TD> <TD align="right"> 85.73 </TD> <TD align="right"> 74.71 </TD> </TR>
  <TR> <TD align="right"> growth_measure </TD> <TD align="right"> 230.24 </TD> <TD align="right"> 194.09 </TD> <TD align="right"> 174.56 </TD> <TD align="right"> 242.25 </TD> <TD align="right"> 219.68 </TD> <TD align="right"> 330.80 </TD> <TD align="right"> 149.25 </TD> <TD align="right"> 293.81 </TD> <TD align="right"> 229.93 </TD> </TR>
  <TR> <TD align="right"> growth_implement </TD> <TD align="right"> 222.82 </TD> <TD align="right"> 440.06 </TD> <TD align="right"> 247.53 </TD> <TD align="right"> 242.68 </TD> <TD align="right"> 228.19 </TD> <TD align="right"> 215.26 </TD> <TD align="right"> 200.92 </TD> <TD align="right"> 267.55 </TD> <TD align="right"> 373.24 </TD> </TR>
  <TR> <TD align="right"> measure_implement </TD> <TD align="right"> 129.78 </TD> <TD align="right"> 104.16 </TD> <TD align="right"> 97.25 </TD> <TD align="right"> 144.06 </TD> <TD align="right"> 96.15 </TD> <TD align="right"> 148.06 </TD> <TD align="right"> 92.06 </TD> <TD align="right"> 163.47 </TD> <TD align="right"> 145.75 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 189.03 </TD> <TD align="right"> 158.18 </TD> <TD align="right"> 152.17 </TD> <TD align="right"> 233.39 </TD> <TD align="right"> 158.29 </TD> <TD align="right"> 209.95 </TD> <TD align="right"> 126.93 </TD> <TD align="right"> 205.58 </TD> <TD align="right"> 199.91 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


