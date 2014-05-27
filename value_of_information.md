






# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   




We consider a Beverton Holt state equation governing population dynamics, \\( f(x,h) = \frac{A x}{1 + B x} \\)



```r
f <- BevHolt
```





With parameters `A` = `1.5` and `B` = `0.05`.



```r
pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
```




Note that the positive stationary root of the model is given by \\( \frac{A-1}{B} \\), or carring capacity `K` = `10`.  
We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  \\( \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\),



```r
price <- 1
c0 <- 0.01
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, 
    c1 = c1)
```




with price `1`, `c0` `0.01` and `c1` `0`. 




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

```
R Version:  R version 2.14.1 (2011-12-22) 

```




**Test with symmetric noise functions**



```r
scenario <- function(policy_g, policy_m, policy_i) {
    
    z_g <- function() rlnorm(1, 0, policy_g)
    z_m <- function() rlnorm(1, 0, policy_m)
    z_i <- function() rlnorm(1, 0, policy_i)
    
    SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, 
        z_g, z_m, z_i, reps = 20000)
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, 
        xT, profit, delta, reward = 0)
}
```






```r
lvl <- 0.2
det <- scenario(0, 0, 0)
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
gmi <- scenario(lvl, lvl, lvl)
```

```
Library ggplot2 loaded.
```






```r
det <- scenario(0.01, 0, 0)
```

```
Library ggplot2 loaded.
```




### plots



```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, det = det$D[, 
    1], g = g$D[, 1], m = m$D[, 1], i = i$D[, 1], gm = gm$D[, 
    1], gi = gi$D[, 1], gmi = gmi$D[, 1]), id = "stock")
ggplot(policy) + geom_point(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7114/7417838132_f5734889ee_o.png) 

```r
ggplot(policy) + geom_smooth(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8005/7417838432_0934048acc_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    g = g$V, m = m$V, i = i$V, gm = gm$V, gi = gi$V, gmi = gmi$V), 
    id = "stock")
ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm6.staticflickr.com/5239/7417838692_dea6f111aa_o.png) 

```r
ggplot(value) + geom_smooth(aes(stock, value, 
    color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7112/7417838928_7c028213c0_o.png) 


## Simulations




```r
simulatereps <- function(opt, true_g, true_m, 
    true_i) {
    z_g <- function() rlnorm(1, 0, true_g)
    z_m <- function() rlnorm(1, 0, true_m)
    z_i <- function() rlnorm(1, 0, true_i)
    
    sims <- lapply(1:100, function(i) {
        ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, 
            opt$D, z_g, z_m, z_i, profit)
    })
    
    # list(sims = sims, opt = opt, true_stochasticity =
    #   c(true_g, true_m, true_i))
    sims
}
```





All cases



```r
policyfn <- list(det = det, g = g, m = m, i = i, 
    gm = gm, gi = gi, gmi = gmi)
noise <- list(s0 = c(0, 0, 0), sg = c(lvl, 0, 
    0), sm = c(0, lvl, 0), si = c(0, 0, lvl), sgm = c(lvl, 
    lvl, 0), sgi = c(lvl, 0, lvl), sgmi = c(lvl, lvl, lvl))
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

![plot of chunk onerep](http://farm6.staticflickr.com/5155/7417843124_84ede34272_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt)
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8011/7417845460_a6c6987f89_o.png) 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, 
    profit)) + facet_wrap(~uncertainty)
```

![The profits made in each time interval of a single replicate, by scenario](http://farm8.staticflickr.com/7108/7417846820_3612cbbc78_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm6.staticflickr.com/5117/7417848506_eb34c83185_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by = uncertainty]
sds <- profits[, sd(V1), by = uncertainty]
```






```r
require(xtable)
uncertainties <- c("deter", "growth", "measure", 
    "implement", "growth+measure", "growth+implement", "all")
print(xtable(matrix(means$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 21 20:33:55 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deter </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deter </TD> <TD align="right"> 33.13 </TD> <TD align="right"> 33.12 </TD> <TD align="right"> 33.13 </TD> <TD align="right"> 32.98 </TD> <TD align="right"> 32.98 </TD> <TD align="right"> 32.93 </TD> <TD align="right"> 32.38 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 35.04 </TD> <TD align="right"> 35.73 </TD> <TD align="right"> 35.15 </TD> <TD align="right"> 35.52 </TD> <TD align="right"> 36.17 </TD> <TD align="right"> 35.65 </TD> <TD align="right"> 34.67 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 30.61 </TD> <TD align="right"> 31.43 </TD> <TD align="right"> 31.66 </TD> <TD align="right"> 30.38 </TD> <TD align="right"> 31.62 </TD> <TD align="right"> 30.78 </TD> <TD align="right"> 30.99 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 32.73 </TD> <TD align="right"> 32.53 </TD> <TD align="right"> 32.66 </TD> <TD align="right"> 32.48 </TD> <TD align="right"> 32.64 </TD> <TD align="right"> 32.33 </TD> <TD align="right"> 32.07 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 32.07 </TD> <TD align="right"> 33.67 </TD> <TD align="right"> 33.84 </TD> <TD align="right"> 32.05 </TD> <TD align="right"> 32.98 </TD> <TD align="right"> 32.83 </TD> <TD align="right"> 33.51 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 34.75 </TD> <TD align="right"> 34.68 </TD> <TD align="right"> 35.16 </TD> <TD align="right"> 33.83 </TD> <TD align="right"> 35.94 </TD> <TD align="right"> 35.70 </TD> <TD align="right"> 34.38 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 31.72 </TD> <TD align="right"> 30.54 </TD> <TD align="right"> 33.01 </TD> <TD align="right"> 31.60 </TD> <TD align="right"> 32.23 </TD> <TD align="right"> 31.57 </TD> <TD align="right"> 32.50 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 21 20:33:55 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deter </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deter </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 5.80 </TD> <TD align="right"> 4.84 </TD> <TD align="right"> 5.57 </TD> <TD align="right"> 6.11 </TD> <TD align="right"> 5.12 </TD> <TD align="right"> 5.45 </TD> <TD align="right"> 6.28 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 3.94 </TD> <TD align="right"> 1.02 </TD> <TD align="right"> 1.01 </TD> <TD align="right"> 4.23 </TD> <TD align="right"> 0.83 </TD> <TD align="right"> 2.19 </TD> <TD align="right"> 3.24 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 0.54 </TD> <TD align="right"> 0.64 </TD> <TD align="right"> 0.60 </TD> <TD align="right"> 0.62 </TD> <TD align="right"> 0.51 </TD> <TD align="right"> 0.64 </TD> <TD align="right"> 0.66 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 7.25 </TD> <TD align="right"> 5.71 </TD> <TD align="right"> 5.35 </TD> <TD align="right"> 7.74 </TD> <TD align="right"> 5.68 </TD> <TD align="right"> 5.80 </TD> <TD align="right"> 5.88 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 5.19 </TD> <TD align="right"> 5.73 </TD> <TD align="right"> 5.78 </TD> <TD align="right"> 6.03 </TD> <TD align="right"> 6.28 </TD> <TD align="right"> 6.30 </TD> <TD align="right"> 5.62 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 6.67 </TD> <TD align="right"> 8.41 </TD> <TD align="right"> 5.87 </TD> <TD align="right"> 6.93 </TD> <TD align="right"> 5.51 </TD> <TD align="right"> 8.04 </TD> <TD align="right"> 7.73 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


