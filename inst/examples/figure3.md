






# Calculating the value of information

 Implements a numerical version of the SDP described in .
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
grid_n <- 200
```


We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `150` on a grid of `200` points, and over an identical discrete grid of possible harvest values.  



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



Do the deterministic exactly,


```r
pdfn <- function(P, s) {
    dunif(P, 1 - s, 1 + s)
}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g = 0.5, pdfn)
det <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward = 0)
```


Determine the policies for each of the scenarios (noise combinations).


```r
lvl <- 0.5
low <- scenario(0.1, 0.1, 0.1)
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


### plots




```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, deterministic = det$D[, 1], all_low = low$D[, 
    1], growth = g$D[, 1], measurement = m$D[, 1], implementation = i$D[, 1]), 
    id = "stock")

ggplot(policy) + geom_point(aes(stock, stock - x_grid[value], color = variable), 
    shape = "+")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8311/7932321418_37f6af0b00_o.png) 

```r
dat <- subset(policy, stock < 140)
dt <- data.table(dat)
linear <- dt[, approx(stock, stock - x_grid[value], xout = seq(1, 150, length = 15)), 
    by = variable]
ggplot(linear) + stat_smooth(aes(x, y, color = variable), degree = 1, se = FALSE, 
    span = 0.3) + xlab("Measured Stock") + ylab("Optimal Expected Escapement")
```

![plot of chunk sethiplots](http://farm9.staticflickr.com/8179/7932321680_2664f358bb_o.png) 


