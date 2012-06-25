



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
function(x, h, p){
  sapply(x, function(x){
         s <- x - h # escapement
         r <- as.numeric(p[1])
         K <- as.numeric(p[2])
         a <- as.numeric(p[3])
         H <- as.numeric(p[4])
         Q <- as.numeric(p[5])
         s * exp(r * (1 - s / K) - a * s ^ (Q - 1) / (s ^ Q + H ^ Q)) 
  })
}
<environment: namespace:pdgControl>
```

```r
curve(0.75 * (1 - x/10), 0, 10)
curve(1 * x^2/(x^3 + 1), 0, 10, add = T, col = "blue")
curve(1.9 * x^2/(x^3 + 1), 0, 10, add = T, col = "red")
```

![plot of chunk showMay](http://farm9.staticflickr.com/8017/7443680468_1dcd08c276_o.png) 



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

![plot of chunk sethiplots](http://farm8.staticflickr.com/7108/7443680954_0e1555a5be_o.png) 

```r

ggplot(policy) + geom_point(aes(stock, x_grid[value], color = variable)) + 
    geom_smooth(aes(stock, x_grid[value], color = variable)) + ylab("harvest")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8023/7443681684_8796c9bb55_o.png) 

```r


value <- melt(data.frame(stock = x_grid, det = det$V, low = low$V, 
    g = g$V, m = m$V, gm = gm$V, gi = gi$V, mi = mi$V, gmi = gmi$V), id = "stock")

ggplot(value) + geom_point(aes(stock, value, color = variable)) + 
    geom_smooth(aes(stock, value, color = variable)) + ylab("Net Present Value")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8156/7443682100_125c7a2d8a_o.png) 


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

```
Error: replacement has length zero
```






```r
sims <- unlist(allcases, recursive = FALSE)
```

```
Error: object 'allcases' not found
```

```r
dat <- melt(sims, id = names(sims[[1]][[1]]))
```

```
Error: object 'sims' not found
```

```r
dt <- data.table(dat)
```

```
Error: object 'dat' not found
```

```r
setnames(dt, c("L2", "L1"), c("reps", "uncertainty"))  # names are nice
```

```
Error: x is not a data.table
```





### Plots 




```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, fishstock)) + 
    geom_line(aes(time, harvest), col = "darkgreen") + facet_wrap(~uncertainty)
```

```
Error: object 'reps' not found
```




This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(subset(dt, fishstock > 0))
```

```
Error: object 'fishstock' not found
```

```r
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + 
    facet_wrap(~uncertainty)
```

```
Error: object 'p1' not found
```







```r
profits <- dt[, sum(profit), by = c("reps", "uncertainty")]
```

```
Error: invalid 'type' (closure) of argument
```

```r
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

```
Error: object 'profits' not found
```




Summary statistics 



```r
means <- profits[, mean(V1), by = uncertainty]
```

```
Error: object 'profits' not found
```

```r
sds <- profits[, sd(V1), by = uncertainty]
```

```
Error: object 'profits' not found
```






```r
require(xtable)
uncertainties <- names(noise)
print(xtable(matrix(means$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

```
Error: object 'means' not found
```

```r
print(xtable(matrix(sds$V1, nrow = length(noise), dimnames = list(uncertainties, 
    uncertainties))), type = "html")
```

```
Error: object 'sds' not found
```





