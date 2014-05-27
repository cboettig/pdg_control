



```
Error: .onLoad failed in loadNamespace() for 'XML', details: call:
dyn.load(file, DLLpath = DLLpath, ...)  error: unable to load shared
object
'/home/cboettig/R/x86_64-redhat-linux-gnu-library/2.15/XML/libs/XML.so':
libxmlsec1.so.1: cannot open shared object file: No such file or directory
```






# Value of information under the May alternative stable states model with uniform noise





Chose the state equation / population dynamics function



```r
f <- May
```





With parameters 



```r
pars <- c(r = 0.75, k = 10, a = 1, H = 1, Q = 3)
K <- 8  # approx
```




Ask R to show us how this function is defined, and plot the transition point using these parameters.



```r
May
```

```
function (x, h, p) 
{
    sapply(x, function(x) {
        s <- max(x - h, 0)
        r <- as.numeric(p[1])
        K <- as.numeric(p[2])
        a <- as.numeric(p[3])
        H <- as.numeric(p[4])
        Q <- as.numeric(p[5])
        s * exp(r * (1 - s/K) - a * s^(Q - 1)/(s^Q + H^Q))
    })
}
<environment: namespace:pdgControl>
```

```r
curve(0.75 * (1 - x/10), 0, 10)
curve(1 * x^2/(x^3 + 1), 0, 10, add = T, col = "blue")
curve(1.9 * x^2/(x^3 + 1), 0, 10, add = T, col = "red")
```

![plot of chunk showMay](figure/showMay.png) 



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




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `12` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




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
R Version:  R version 2.15.0 (2012-03-30) 

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

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, stock - x_grid[value], color = variable)) + ylab("escapement")
```

![plot of chunk sethiplots](figure/sethiplots1.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, x_grid[value], color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](figure/sethiplots2.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](figure/sethiplots3.png) 


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
noise <- list(det = c(0, 0, 0), low = c(0.1, 0.1, 0.1), grow = c(lvl, 
    0, 0), meas = c(0, lvl, 0), imp = c(0, 0, lvl), gro_meas = c(lvl, lvl, 0), 
    gro_imp = c(lvl, 0, lvl), meas_imp = c(0, lvl, lvl), all = c(lvl, lvl, lvl))
```







```r
allcases <- lapply(policyfn, function(policyfn_i) {
    lapply(noise, function(noise_i) {
        try(simulatereps(policyfn_i, noise_i[1], noise_i[2], noise_i[3]))
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

![plot of chunk onerep](figure/onerep.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](figure/stock.png) 





```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](figure/profits.png) 


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

<!-- html table generated in R 2.15.0 by xtable 1.7-0 package -->
<!-- Tue Jun 26 22:45:26 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> grow </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 27.74 </TD> <TD align="right"> 30.18 </TD> <TD align="right"> 29.82 </TD> <TD align="right"> 29.17 </TD> <TD align="right"> 29.67 </TD> <TD align="right"> 28.94 </TD> <TD align="right"> 30.18 </TD> <TD align="right"> 28.73 </TD> <TD align="right"> 30.18 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 27.23 </TD> <TD align="right"> 29.56 </TD> <TD align="right"> 29.59 </TD> <TD align="right"> 29.44 </TD> <TD align="right"> 29.11 </TD> <TD align="right"> 29.01 </TD> <TD align="right"> 29.70 </TD> <TD align="right"> 28.99 </TD> <TD align="right"> 29.64 </TD> </TR>
  <TR> <TD align="right"> grow </TD> <TD align="right"> 26.23 </TD> <TD align="right"> 28.78 </TD> <TD align="right"> 27.78 </TD> <TD align="right"> 26.65 </TD> <TD align="right"> 27.66 </TD> <TD align="right"> 27.88 </TD> <TD align="right"> 27.64 </TD> <TD align="right"> 28.83 </TD> <TD align="right"> 27.62 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 26.23 </TD> <TD align="right"> 24.79 </TD> <TD align="right"> 26.10 </TD> <TD align="right"> 26.25 </TD> <TD align="right"> 25.60 </TD> <TD align="right"> 26.40 </TD> <TD align="right"> 24.81 </TD> <TD align="right"> 26.62 </TD> <TD align="right"> 26.19 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 26.96 </TD> <TD align="right"> 29.27 </TD> <TD align="right"> 28.86 </TD> <TD align="right"> 29.03 </TD> <TD align="right"> 28.73 </TD> <TD align="right"> 28.77 </TD> <TD align="right"> 29.26 </TD> <TD align="right"> 28.74 </TD> <TD align="right"> 29.20 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 23.81 </TD> <TD align="right"> 21.62 </TD> <TD align="right"> 23.45 </TD> <TD align="right"> 24.37 </TD> <TD align="right"> 22.94 </TD> <TD align="right"> 23.00 </TD> <TD align="right"> 21.85 </TD> <TD align="right"> 24.07 </TD> <TD align="right"> 23.23 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 27.68 </TD> <TD align="right"> 26.34 </TD> <TD align="right"> 28.50 </TD> <TD align="right"> 26.71 </TD> <TD align="right"> 27.69 </TD> <TD align="right"> 26.62 </TD> <TD align="right"> 26.94 </TD> <TD align="right"> 27.09 </TD> <TD align="right"> 26.57 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 24.34 </TD> <TD align="right"> 22.55 </TD> <TD align="right"> 23.53 </TD> <TD align="right"> 24.40 </TD> <TD align="right"> 23.74 </TD> <TD align="right"> 25.11 </TD> <TD align="right"> 22.54 </TD> <TD align="right"> 24.93 </TD> <TD align="right"> 24.11 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 21.95 </TD> <TD align="right"> 20.13 </TD> <TD align="right"> 20.57 </TD> <TD align="right"> 21.17 </TD> <TD align="right"> 20.12 </TD> <TD align="right"> 22.27 </TD> <TD align="right"> 19.99 </TD> <TD align="right"> 22.61 </TD> <TD align="right"> 21.09 </TD> </TR>
   </TABLE>


```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

<!-- html table generated in R 2.15.0 by xtable 1.7-0 package -->
<!-- Tue Jun 26 22:45:26 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> det </TH> <TH> low </TH> <TH> grow </TH> <TH> meas </TH> <TH> imp </TH> <TH> gro_meas </TH> <TH> gro_imp </TH> <TH> meas_imp </TH> <TH> all </TH>  </TR>
  <TR> <TD align="right"> det </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> low </TD> <TD align="right"> 1.91 </TD> <TD align="right"> 1.58 </TD> <TD align="right"> 1.73 </TD> <TD align="right"> 1.66 </TD> <TD align="right"> 1.84 </TD> <TD align="right"> 1.60 </TD> <TD align="right"> 1.73 </TD> <TD align="right"> 1.77 </TD> <TD align="right"> 1.45 </TD> </TR>
  <TR> <TD align="right"> grow </TD> <TD align="right"> 8.58 </TD> <TD align="right"> 9.09 </TD> <TD align="right"> 9.19 </TD> <TD align="right"> 9.28 </TD> <TD align="right"> 9.03 </TD> <TD align="right"> 7.72 </TD> <TD align="right"> 9.24 </TD> <TD align="right"> 8.68 </TD> <TD align="right"> 8.87 </TD> </TR>
  <TR> <TD align="right"> meas </TD> <TD align="right"> 1.55 </TD> <TD align="right"> 1.99 </TD> <TD align="right"> 1.35 </TD> <TD align="right"> 1.18 </TD> <TD align="right"> 1.63 </TD> <TD align="right"> 1.23 </TD> <TD align="right"> 1.86 </TD> <TD align="right"> 1.16 </TD> <TD align="right"> 1.23 </TD> </TR>
  <TR> <TD align="right"> imp </TD> <TD align="right"> 0.84 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 1.03 </TD> <TD align="right"> 0.35 </TD> <TD align="right"> 0.97 </TD> <TD align="right"> 0.64 </TD> <TD align="right"> 0.94 </TD> <TD align="right"> 0.45 </TD> <TD align="right"> 0.87 </TD> </TR>
  <TR> <TD align="right"> gro_meas </TD> <TD align="right"> 7.37 </TD> <TD align="right"> 6.80 </TD> <TD align="right"> 7.82 </TD> <TD align="right"> 6.73 </TD> <TD align="right"> 7.47 </TD> <TD align="right"> 7.69 </TD> <TD align="right"> 6.93 </TD> <TD align="right"> 7.54 </TD> <TD align="right"> 7.38 </TD> </TR>
  <TR> <TD align="right"> gro_imp </TD> <TD align="right"> 8.96 </TD> <TD align="right"> 8.31 </TD> <TD align="right"> 7.64 </TD> <TD align="right"> 8.29 </TD> <TD align="right"> 7.76 </TD> <TD align="right"> 7.94 </TD> <TD align="right"> 8.47 </TD> <TD align="right"> 8.98 </TD> <TD align="right"> 7.59 </TD> </TR>
  <TR> <TD align="right"> meas_imp </TD> <TD align="right"> 2.97 </TD> <TD align="right"> 4.41 </TD> <TD align="right"> 3.99 </TD> <TD align="right"> 2.32 </TD> <TD align="right"> 4.28 </TD> <TD align="right"> 2.45 </TD> <TD align="right"> 4.26 </TD> <TD align="right"> 2.91 </TD> <TD align="right"> 3.30 </TD> </TR>
  <TR> <TD align="right"> all </TD> <TD align="right"> 8.12 </TD> <TD align="right"> 7.04 </TD> <TD align="right"> 6.80 </TD> <TD align="right"> 7.41 </TD> <TD align="right"> 6.93 </TD> <TD align="right"> 7.79 </TD> <TD align="right"> 7.50 </TD> <TD align="right"> 7.41 </TD> <TD align="right"> 6.45 </TD> </TR>
   </TABLE>




