






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
    
    z_g <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_g
    z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_m
    z_i <- function() 1 + (2 * runif(1, 0, 1) - 1) * true_i
    
    sims <- lapply(1:100, function(i) {
        ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, 
            opt$D, z_g, z_m, z_i, profit)
    })
    
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

![plot of chunk onerep](http://farm8.staticflickr.com/7262/7416382244_b42ea9aae1_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt)
p1 + geom_line(aes(time, fishstock, group = reps), 
    alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8003/7416386694_7a2beb2e6f_o.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7125/7416389440_b7bdf89b4f_o.png) 


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
<!-- Thu Jun 21 14:34:05 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deterministic </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deterministic </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> <TD align="right"> 674.24 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 673.86 </TD> <TD align="right"> 672.03 </TD> <TD align="right"> 674.80 </TD> <TD align="right"> 672.41 </TD> <TD align="right"> 675.34 </TD> <TD align="right"> 673.48 </TD> <TD align="right"> 674.64 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 668.34 </TD> <TD align="right"> 668.35 </TD> <TD align="right"> 668.49 </TD> <TD align="right"> 668.40 </TD> <TD align="right"> 668.90 </TD> <TD align="right"> 668.50 </TD> <TD align="right"> 668.48 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 672.65 </TD> <TD align="right"> 672.22 </TD> <TD align="right"> 672.62 </TD> <TD align="right"> 671.83 </TD> <TD align="right"> 672.40 </TD> <TD align="right"> 672.56 </TD> <TD align="right"> 672.22 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 666.84 </TD> <TD align="right"> 672.25 </TD> <TD align="right"> 670.16 </TD> <TD align="right"> 667.48 </TD> <TD align="right"> 667.78 </TD> <TD align="right"> 668.48 </TD> <TD align="right"> 672.21 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 673.36 </TD> <TD align="right"> 668.46 </TD> <TD align="right"> 670.25 </TD> <TD align="right"> 671.88 </TD> <TD align="right"> 672.58 </TD> <TD align="right"> 671.20 </TD> <TD align="right"> 673.67 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 665.96 </TD> <TD align="right"> 667.67 </TD> <TD align="right"> 664.82 </TD> <TD align="right"> 665.86 </TD> <TD align="right"> 667.26 </TD> <TD align="right"> 667.77 </TD> <TD align="right"> 666.34 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = 7, dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.14.1 by xtable 1.7-0 package -->
<!-- Thu Jun 21 14:34:05 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> deterministic </TH> <TH> growth </TH> <TH> measure </TH> <TH> implement </TH> <TH> growth+measure </TH> <TH> growth+implement </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> deterministic </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> growth </TD> <TD align="right"> 22.64 </TD> <TD align="right"> 20.26 </TD> <TD align="right"> 20.05 </TD> <TD align="right"> 18.94 </TD> <TD align="right"> 19.55 </TD> <TD align="right"> 20.95 </TD> <TD align="right"> 20.65 </TD> </TR>
  <TR> <TD align="right"> measure </TD> <TD align="right"> 2.64 </TD> <TD align="right"> 2.65 </TD> <TD align="right"> 2.50 </TD> <TD align="right"> 2.66 </TD> <TD align="right"> 2.34 </TD> <TD align="right"> 2.50 </TD> <TD align="right"> 2.35 </TD> </TR>
  <TR> <TD align="right"> implement </TD> <TD align="right"> 2.49 </TD> <TD align="right"> 2.29 </TD> <TD align="right"> 2.44 </TD> <TD align="right"> 2.73 </TD> <TD align="right"> 2.43 </TD> <TD align="right"> 2.48 </TD> <TD align="right"> 2.21 </TD> </TR>
  <TR> <TD align="right"> growth+measure </TD> <TD align="right"> 22.89 </TD> <TD align="right"> 22.07 </TD> <TD align="right"> 22.25 </TD> <TD align="right"> 21.92 </TD> <TD align="right"> 21.28 </TD> <TD align="right"> 19.32 </TD> <TD align="right"> 21.11 </TD> </TR>
  <TR> <TD align="right"> growth+implement </TD> <TD align="right"> 18.83 </TD> <TD align="right"> 22.31 </TD> <TD align="right"> 20.58 </TD> <TD align="right"> 22.80 </TD> <TD align="right"> 21.05 </TD> <TD align="right"> 21.30 </TD> <TD align="right"> 20.11 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 23.53 </TD> <TD align="right"> 21.18 </TD> <TD align="right"> 20.88 </TD> <TD align="right"> 21.54 </TD> <TD align="right"> 20.23 </TD> <TD align="right"> 20.53 </TD> <TD align="right"> 19.52 </TD> </TR>
   </TABLE>





# References

<p>Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
&ldquo;Fishery Management Under Multiple Uncertainty.&rdquo;
<EM>Journal of Environmental Economics And Management</EM>, <B>50</B>.
ISSN 00950696, <a href="http://dx.doi.org/10.1016/j.jeem.2004.11.005">http://dx.doi.org/10.1016/j.jeem.2004.11.005</a>.


