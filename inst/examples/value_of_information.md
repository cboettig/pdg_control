






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
det <- scenario(0, 0, 0)
```

```
Library ggplot2 loaded.
```

```r
g <- scenario(0.2, 0, 0)
```

```
Library ggplot2 loaded.
```

```r
m <- scenario(0, 0.2, 0)
```

```
Library ggplot2 loaded.
```

```r
i <- scenario(0, 0, 0.2)
```

```
Library ggplot2 loaded.
```

```r
gm <- scenario(0.2, 0.2, 0)
```

```
Library ggplot2 loaded.
```

```r
gi <- scenario(0.2, 0, 0.2)
```

```
Library ggplot2 loaded.
```

```r
gmi <- scenario(0.2, 0.2, 0.2)
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
    1], g = g$D[, 1], m = m$D[, 1], gm = gm$D[, 1], gi = gi$D[, 
    1], gmi = gmi$D[, 1]), id = "stock")
ggplot(policy) + geom_point(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7095/7395320626_34f2212554_o.png) 

```r
ggplot(policy) + geom_smooth(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7216/7395321048_556df76031_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, gmi = gmi$V), 
    id = "stock")
ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8006/7395321480_7613246778_o.png) 

```r
ggplot(value) + geom_smooth(aes(stock, value, 
    color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8141/7395322218_7212067722_o.png) 


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
noise <- list(s0 = c(0, 0, 0), sg = c(0.2, 0, 
    0), sm = c(0, 0.2, 0), si = c(0, 0, 0.2), sgm = c(0.2, 
    0.2, 0), sgi = c(0.2, 0, 0.2), sgmi = c(0.2, 0.2, 0.2))
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

![plot of chunk onerep](http://farm6.staticflickr.com/5333/7395548324_acedc3251d_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt)
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm6.staticflickr.com/5275/7395552194_95bbab1150_o.png) 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, 
    profit)) + facet_wrap(~uncertainty)
```

![The profits made in each time interval of a single replicate, by scenario](http://farm8.staticflickr.com/7241/7395628756_6dd34b9d28_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm9.staticflickr.com/8168/7395556128_1e80a34153_o.png) 


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
<!-- Mon Jun 18 11:22:39 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deter </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deter </TD> <TD align="right"> 33.13 </TD> <TD align="right"> 33.13 </TD> <TD align="right"> 33.13 </TD> <TD align="right"> 33.07 </TD> <TD align="right"> 32.98 </TD> <TD align="right"> 32.62 </TD> <TD align="right"> 32.91 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 36.02 </TD> <TD align="right"> 35.41 </TD> <TD align="right"> 33.83 </TD> <TD align="right"> 36.73 </TD> <TD align="right"> 35.46 </TD> <TD align="right"> 35.43 </TD> <TD align="right"> 35.86 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 30.72 </TD> <TD align="right"> 30.25 </TD> <TD align="right"> 31.37 </TD> <TD align="right"> 30.84 </TD> <TD align="right"> 31.64 </TD> <TD align="right"> 30.97 </TD> <TD align="right"> 31.44 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 32.66 </TD> <TD align="right"> 32.51 </TD> <TD align="right"> 32.65 </TD> <TD align="right"> 32.21 </TD> <TD align="right"> 32.64 </TD> <TD align="right"> 31.97 </TD> <TD align="right"> 31.95 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 31.55 </TD> <TD align="right"> 32.23 </TD> <TD align="right"> 34.26 </TD> <TD align="right"> 32.96 </TD> <TD align="right"> 32.60 </TD> <TD align="right"> 34.06 </TD> <TD align="right"> 33.77 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 34.51 </TD> <TD align="right"> 34.73 </TD> <TD align="right"> 35.12 </TD> <TD align="right"> 35.38 </TD> <TD align="right"> 35.57 </TD> <TD align="right"> 35.47 </TD> <TD align="right"> 35.64 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 31.19 </TD> <TD align="right"> 31.02 </TD> <TD align="right"> 31.71 </TD> <TD align="right"> 29.65 </TD> <TD align="right"> 32.47 </TD> <TD align="right"> 32.58 </TD> <TD align="right"> 32.23 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Mon Jun 18 11:22:39 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deter </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deter </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 5.70 </TD> <TD align="right"> 5.52 </TD> <TD align="right"> 5.61 </TD> <TD align="right"> 5.95 </TD> <TD align="right"> 7.11 </TD> <TD align="right"> 5.63 </TD> <TD align="right"> 5.36 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 3.14 </TD> <TD align="right"> 4.61 </TD> <TD align="right"> 2.04 </TD> <TD align="right"> 3.12 </TD> <TD align="right"> 1.45 </TD> <TD align="right"> 2.64 </TD> <TD align="right"> 0.89 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 0.58 </TD> <TD align="right"> 0.73 </TD> <TD align="right"> 0.57 </TD> <TD align="right"> 2.34 </TD> <TD align="right"> 0.48 </TD> <TD align="right"> 2.30 </TD> <TD align="right"> 2.37 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 7.64 </TD> <TD align="right"> 5.78 </TD> <TD align="right"> 4.88 </TD> <TD align="right"> 6.74 </TD> <TD align="right"> 5.55 </TD> <TD align="right"> 5.46 </TD> <TD align="right"> 5.74 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 5.30 </TD> <TD align="right"> 5.55 </TD> <TD align="right"> 5.18 </TD> <TD align="right"> 5.50 </TD> <TD align="right"> 6.41 </TD> <TD align="right"> 6.33 </TD> <TD align="right"> 5.00 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 7.51 </TD> <TD align="right"> 7.22 </TD> <TD align="right"> 6.62 </TD> <TD align="right"> 8.60 </TD> <TD align="right"> 6.86 </TD> <TD align="right"> 5.51 </TD> <TD align="right"> 6.96 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


