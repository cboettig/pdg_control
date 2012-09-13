




## Generate Figure 3 of Sethi et al. 2006




Chose the state equation / population dynamics function


```r
f <- function(x, h, p) {
    sapply(x, function(x) {
        S = max(x - h, 0)
        p[1] * S * (1 - S/p[2]) + S
    })
}
```


With parameters 


```r
r <- 1
K <- 100
pars <- c(r, K)
```


`r` = `1` and `K` = `100`.

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
sigma_g <- 0.5
sigma_m <- 0.5
sigma_i <- 0.5
```


We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  

# Scenarios: 


```r
require(cubature)
# Confirm that these give the same value, and time performance
system.time(a <- sapply(x_grid, function(x) int_f(f, x, 1, 0.1, 0.1, pars)))
```

```
   user  system elapsed 
  3.562   0.003   3.576 
```

```r
system.time(b <- sapply(x_grid, function(x) F(x, 1, 0.1, 0.1, pars)))
```

```
   user  system elapsed 
  0.006   0.000   0.007 
```



```r
sdp <- SDP_uniform(f, pars, x_grid, h_grid, sigma_g = sigma_g, pdfn = function(P, 
    s) dunif(P, 1 - s, 1 + s), sigma_m = sigma_m, sigma_i = sigma_i, F)
```



```r
uniform <- find_dp_optim(sdp, x_grid, h_grid, OptTime, xT, profit, delta, reward = 0)
```


Do the deterministic exactly,


```r
pdfn <- function(P, s) {
    dunif(P, 1 - s, 1 + s)
}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g = sigma_g, 
    pdfn)
det <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward = 0)
```





Determine the policies for each of the scenarios (noise combinations).




### plots




```r
require(reshape2)
# policy <- melt(data.frame(stock = x_grid, deterministic=det$D[,1],
# uniform = uniform$D[,1], mc=opt$D[,1]), id = 'stock')
policy <- melt(data.frame(stock = x_grid, deterministic = det$D[, 1], uniform = uniform$D[, 
    1]), id = "stock")

ggplot(policy) + geom_jitter(aes(stock, stock - x_grid[value], color = variable), 
    shape = "+")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8176/7981119045_8dcd1320b5_o.png) 

```r
dat <- subset(policy, stock < 140)
dt <- data.table(dat)
linear <- dt[, approx(stock, stock - x_grid[value], xout = seq(1, 150, length = 15)), 
    by = variable]
ggplot(linear) + stat_smooth(aes(x, y, color = variable), degree = 1, se = FALSE, 
    span = 0.3) + xlab("Measured Stock") + ylab("Optimal Expected Escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8445/7981121116_352840bbde_o.png) 


