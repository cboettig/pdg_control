---
layout: page
---









We consider adding additional sources of noise through observation error and implementation error to the orginal model of Reed (1979), which considers only growth error.  This follows in the line of exploration taken by Clark and Kirkwood (1987), Roughgarden (1996) and Sethi (2006). 

In adding these two sources of error, our state variable becomes our *measured* (or believed) stock $x$, rather than the true stock size $y$, and our control variable becomes the quota $q$ set, rather than the realized harvest $h$.  

## Analytical integrals for Sethi (2006) calcuations

The central calculation to accomodating additional sources of uncertainty is the determination of the transition matrix; the probability of going from state $x_t$ at time $t$ to state $x_{t+1}$ in the next interval, for each possible value of the control variable (quota).  We discritize the statespace onto a fine grid of values $x$ and $h$ to seek a stochastic dynamic programming solution.  



```r
xmin <- 0
xmax <- 150
grid_n <- 150
```


We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `150` on a grid of `150` points, and over an identical discrete grid of possible harvest values.  



```r
x_grid <- seq(xmin, xmax, length = grid_n)
h_grid <- x_grid
```



With purely growth noise, a row of the matrix is given by the probability density of the growth noise with mean $f(x,h)$, that is, the probability of tranistion from $x$ to $x_{t+1}$ at harvest $h$ is given by the function $P(x_t+1, f(x,h))$ for some probability density function $P$. 



We assume a simple function like logistic growth, \\(f(s) = r s (1-s/K) + s \\), where \\(s = x - h\\),


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


## Numerical implementation


We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  $\Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x}$, conditioned on $h > x$ and $x > 0$


```r
price <- 1
c0 <- 0
c1 <- 0
profit <- profit_harvest(price = price, c0 = c0, c1 = c1)
```


with price = `1`, `c0` = `0` and `c1` = `0`. 


Additional parameters


```r
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- 0.5
sigma_m <- 0.5
sigma_i <- 0.5
```



The uniform distribution must be careful for noise sizes smaller than binwidth


```r
pdfn = function(P, mu, s) {
    if (mu == 0) 
        if (P == 0) 
            1 else 0 else dunif(P, mu * (1 - s), mu * (1 + s))
}

ln <- function(P, mu, s) dlnorm(P, log(mu), s)
```


We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  

# Scenarios: 



```r
sdp <- SDP_multiple_uncertainty(f, pars, x_grid, h_grid, sigma_g, 
    pdfn, sigma_m, sigma_i)
```



```r
mult <- find_dp_optim(sdp, x_grid, h_grid, OptTime, xT, profit, delta, 
    reward = 0)
```


Compare to the former method:


```r
det_mat <- SDP_multiple_uncertainty(f, pars, x_grid, h_grid, 0, pdfn, 
    0, 0)
det <- find_dp_optim(det_mat, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward = 0)
```




### plots




```r
require(reshape2)

policy <- melt(data.frame(stock = x_grid, deterministic = det$D[, 
    1], new = mult$D[, 1]), id = "stock")
ggplot(subset(policy, stock < 120)) + geom_jitter(aes(stock, stock - 
    x_grid[value], color = variable), shape = "+")
```

![plot of chunk sethiplots](figure/sethiplots1.png) 

```r

dat <- subset(policy, stock < 120)
dt <- data.table(dat)
linear <- dt[, approx(stock, stock - x_grid[value], xout = seq(1, 
    120, length = 15)), by = variable]
ggplot(linear) + stat_smooth(aes(x, y, color = variable), degree = 1, 
    se = FALSE, span = 0.3) + xlab("Measured Stock") + ylab("Optimal Expected Escapement")
```

![plot of chunk sethiplots](figure/sethiplots2.png) 


