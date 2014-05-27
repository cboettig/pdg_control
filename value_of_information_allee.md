






# Calculating the value of information with highly nonlinear value function

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty, but under strong allee dynamics 




We consider the state dynamics 



```r
f <- RickerAllee
K <- 4
xT <- 2  # final value, also allee threshold
pars <- c(1.5, K, xT)
f
```

```
function (x, h, p) 
{
    sapply(x, function(x) {
        x <- max(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])
    })
}
<environment: namespace:pdgControl>
```




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




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `6` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




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
g <- scenario(0.1, 0, 0)
```

```
Library ggplot2 loaded.
```

```r
m <- scenario(0, 0.1, 0)
```

```
Library ggplot2 loaded.
```

```r
i <- scenario(0, 0, 0.1)
```

```
Library ggplot2 loaded.
```

```r
gm <- scenario(0.1, 0.1, 0)
```

```
Library ggplot2 loaded.
```

```r
gi <- scenario(0.1, 0, 0.1)
```

```
Library ggplot2 loaded.
```

```r
gmi <- scenario(0.1, 0.1, 0.1)
```

```
Library ggplot2 loaded.
```






```r
det <- scenario(0.001, 0, 0)
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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7075/7395812422_7104273350_o.png) 

```r
ggplot(policy) + geom_smooth(aes(stock, stock - 
    x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8026/7395813164_3d321c9d2a_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, gmi = gmi$V), 
    id = "stock")
ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7242/7395813868_e2c7f6c266_o.png) 

```r
ggplot(value) + geom_smooth(aes(stock, value, 
    color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm8.staticflickr.com/7220/7395814590_53afcce0f5_o.png) 


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

![plot of chunk onerep](http://farm6.staticflickr.com/5463/7395822046_7f264bdc3d_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt)
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7217/7395826070_02f12b0684_o.png) 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, 
    profit)) + facet_wrap(~uncertainty)
```

![The profits made in each time interval of a single replicate, by scenario](http://farm6.staticflickr.com/5198/7395827978_cbae7381c5_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7092/7395830152_b25eb953db_o.png) 


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
    uncertainties)), caption = "Mean realized net present value over simulations"), 
    type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Mon Jun 18 11:30:07 2012 -->
<TABLE border=1>
<CAPTION ALIGN="bottom"> Mean realized net present value over simulations </CAPTION>
<TR> <TH>  </TH> <TH> deter </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deter </TD> <TD align="right"> 10.92 </TD> <TD align="right"> 10.91 </TD> <TD align="right"> 10.92 </TD> <TD align="right"> 10.86 </TD> <TD align="right"> 10.88 </TD> <TD align="right"> 10.64 </TD> <TD align="right"> 10.68 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 9.73 </TD> <TD align="right"> 9.43 </TD> <TD align="right"> 9.94 </TD> <TD align="right"> 10.68 </TD> <TD align="right"> 9.69 </TD> <TD align="right"> 10.47 </TD> <TD align="right"> 11.02 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 5.41 </TD> <TD align="right"> 6.10 </TD> <TD align="right"> 5.81 </TD> <TD align="right"> 5.75 </TD> <TD align="right"> 5.84 </TD> <TD align="right"> 6.33 </TD> <TD align="right"> 6.52 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 10.63 </TD> <TD align="right"> 10.62 </TD> <TD align="right"> 10.62 </TD> <TD align="right"> 10.51 </TD> <TD align="right"> 10.55 </TD> <TD align="right"> 10.37 </TD> <TD align="right"> 10.36 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 5.27 </TD> <TD align="right"> 5.64 </TD> <TD align="right"> 4.68 </TD> <TD align="right"> 5.73 </TD> <TD align="right"> 6.07 </TD> <TD align="right"> 5.98 </TD> <TD align="right"> 6.76 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 9.53 </TD> <TD align="right"> 9.33 </TD> <TD align="right"> 8.83 </TD> <TD align="right"> 9.64 </TD> <TD align="right"> 9.23 </TD> <TD align="right"> 10.37 </TD> <TD align="right"> 10.47 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 4.94 </TD> <TD align="right"> 5.63 </TD> <TD align="right"> 5.39 </TD> <TD align="right"> 5.29 </TD> <TD align="right"> 5.98 </TD> <TD align="right"> 6.25 </TD> <TD align="right"> 5.09 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties)), caption = "Standard deviation in realized net present value over simulations"), 
    type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Mon Jun 18 11:30:07 2012 -->
<TABLE border=1>
<CAPTION ALIGN="bottom"> Standard deviation in realized net present value over simulations </CAPTION>
<TR> <TH>  </TH> <TH> deter </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deter </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 5.08 </TD> <TD align="right"> 4.06 </TD> <TD align="right"> 4.46 </TD> <TD align="right"> 4.04 </TD> <TD align="right"> 4.83 </TD> <TD align="right"> 4.58 </TD> <TD align="right"> 4.69 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 1.52 </TD> <TD align="right"> 1.95 </TD> <TD align="right"> 1.84 </TD> <TD align="right"> 1.83 </TD> <TD align="right"> 1.72 </TD> <TD align="right"> 2.06 </TD> <TD align="right"> 2.04 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 0.35 </TD> <TD align="right"> 0.33 </TD> <TD align="right"> 0.36 </TD> <TD align="right"> 0.39 </TD> <TD align="right"> 0.36 </TD> <TD align="right"> 0.37 </TD> <TD align="right"> 0.34 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 2.30 </TD> <TD align="right"> 2.87 </TD> <TD align="right"> 1.80 </TD> <TD align="right"> 3.00 </TD> <TD align="right"> 2.86 </TD> <TD align="right"> 2.86 </TD> <TD align="right"> 3.40 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 4.45 </TD> <TD align="right"> 4.33 </TD> <TD align="right"> 4.26 </TD> <TD align="right"> 4.41 </TD> <TD align="right"> 4.33 </TD> <TD align="right"> 4.39 </TD> <TD align="right"> 4.41 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 1.94 </TD> <TD align="right"> 2.84 </TD> <TD align="right"> 2.32 </TD> <TD align="right"> 2.31 </TD> <TD align="right"> 3.04 </TD> <TD align="right"> 3.35 </TD> <TD align="right"> 2.46 </TD> </TR>
   </TABLE>







# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


