






# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   




Chose the state equation / population dynamics function



```r
f <- function(x, h, p) {
    S = x - h
    p[1] * S * (1 - S/p[2]) + S
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






```r
lvl <- 0.1
det <- scenario(0.01, 0, 0)
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
im <- scenario(0, lvl, lvl)
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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7259/7416278218_b7b4271ce3_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], 
    color = variable)) + geom_smooth(aes(stock, x_grid[value], 
    color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8164/7416278820_8704828f66_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, gmi = gmi$V), 
    id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8005/7416279222_c4e163b8c7_o.png) 


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

![plot of chunk onerep](http://farm8.staticflickr.com/7266/7415825718_b52a7acc4d_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt)
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm6.staticflickr.com/5197/7415831018_e4f7b36dd6_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm6.staticflickr.com/5075/7415834216_c6bce3e201_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by = uncertainty]
sds <- profits[, sd(V1), by = uncertainty]
```






```r
require(xtable)
uncertainties <- c("deterministic", "growth", 
    "measure", "implement", "growth+measure", "growth+implement", 
    "all")
print(xtable(matrix(means$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 21 14:11:45 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deterministic </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deterministic </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 684.77 </TD> <TD align="right"> 687.36 </TD> <TD align="right"> 680.98 </TD> <TD align="right"> 688.88 </TD> <TD align="right"> 686.54 </TD> <TD align="right"> 687.39 </TD> <TD align="right"> 677.78 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 659.35 </TD> <TD align="right"> 659.19 </TD> <TD align="right"> 657.10 </TD> <TD align="right"> 657.44 </TD> <TD align="right"> 658.29 </TD> <TD align="right"> 658.09 </TD> <TD align="right"> 658.41 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 670.45 </TD> <TD align="right"> 670.92 </TD> <TD align="right"> 670.80 </TD> <TD align="right"> 670.37 </TD> <TD align="right"> 670.24 </TD> <TD align="right"> 670.34 </TD> <TD align="right"> 670.30 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 670.09 </TD> <TD align="right"> 669.03 </TD> <TD align="right"> 662.09 </TD> <TD align="right"> 663.28 </TD> <TD align="right"> 665.16 </TD> <TD align="right"> 664.49 </TD> <TD align="right"> 664.75 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 676.98 </TD> <TD align="right"> 675.74 </TD> <TD align="right"> 680.50 </TD> <TD align="right"> 680.06 </TD> <TD align="right"> 678.51 </TD> <TD align="right"> 675.86 </TD> <TD align="right"> 678.70 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 666.08 </TD> <TD align="right"> 664.86 </TD> <TD align="right"> 666.87 </TD> <TD align="right"> 659.18 </TD> <TD align="right"> 659.58 </TD> <TD align="right"> 662.80 </TD> <TD align="right"> 658.50 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 21 14:11:45 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deterministic </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deterministic </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 36.26 </TD> <TD align="right"> 36.28 </TD> <TD align="right"> 38.06 </TD> <TD align="right"> 37.42 </TD> <TD align="right"> 32.84 </TD> <TD align="right"> 34.16 </TD> <TD align="right"> 42.24 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 5.64 </TD> <TD align="right"> 5.49 </TD> <TD align="right"> 6.75 </TD> <TD align="right"> 7.33 </TD> <TD align="right"> 6.15 </TD> <TD align="right"> 5.61 </TD> <TD align="right"> 5.56 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 3.64 </TD> <TD align="right"> 3.40 </TD> <TD align="right"> 3.66 </TD> <TD align="right"> 4.10 </TD> <TD align="right"> 4.60 </TD> <TD align="right"> 3.92 </TD> <TD align="right"> 3.34 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 37.59 </TD> <TD align="right"> 38.37 </TD> <TD align="right"> 29.85 </TD> <TD align="right"> 38.17 </TD> <TD align="right"> 35.71 </TD> <TD align="right"> 35.61 </TD> <TD align="right"> 38.08 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 38.81 </TD> <TD align="right"> 35.78 </TD> <TD align="right"> 38.83 </TD> <TD align="right"> 38.85 </TD> <TD align="right"> 41.09 </TD> <TD align="right"> 36.06 </TD> <TD align="right"> 34.98 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 38.60 </TD> <TD align="right"> 37.18 </TD> <TD align="right"> 32.69 </TD> <TD align="right"> 40.24 </TD> <TD align="right"> 34.89 </TD> <TD align="right"> 38.93 </TD> <TD align="right"> 37.73 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


