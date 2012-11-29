










```r
seed <- 123  # Random seed (replicable results)
delta <- 0.05  # economic discounting rate
OptTime <- 100  # stopping time
gridsize <- 50  # grid size for fish stock and harvest rate (discretized population)
sigma_g <- 0.2  # Noise in population growth
reward <- 0  # bonus for satisfying the boundary condition
z_g <- function() rlnorm(1, 0, sigma_g)  # mean 1
z_m <- function() 1  # No measurement noise,
z_i <- function() 1  # No implemenation noise
f <- BevHolt  # Select the state equation
pars <- c(1.5, 0.05)  # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference
xT <- 0  # boundary conditions
x0 <- K
profit <- profit_harvest(price = 10, c0 = 0, c1 = 0)
x_grid <- seq(0.01, 1.2 * K, length = gridsize)
h_grid <- seq(0.01, 0.8 * K, length = gridsize)
```




```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward = reward)
```





```r
L1 <- function(c2) function(h, h_prev) c2 * abs(h - h_prev)
free_increase <- function(c2) function(h, h_prev) c2 * abs(min(h - h_prev, 0))  # increasing harvest is free
free_decrease <- function(c2) function(h, h_prev) c2 * max(h - h_prev, 0)  # decreasing harvest is free
fixed <- function(c2) function(h, h_prev) c2 * as.numeric(!(h == h_prev))
L2 <- function(c2) function(h, h_prev) c2 * (h - h_prev)^2
none <- function(h, h_prev) 0
penaltyfns <- list(L2 = L2, L1 = L1, free_decrease = free_decrease, fixed = fixed, 
    free_increase = free_increase)
c2 <- seq(0, 10, length.out = 30)
```


## Apples to Apples levels

This can take a while, so we use explicit parallelization, 

```r
require(snowfall)
sfInit(cpu = 8, parallel = T)
```

```
R Version:  R version 2.15.2 (2012-10-26) 

```

```r
sfLibrary(pdgControl)
```

```
Library pdgControl loaded.
```

```r
sfExportAll()
```


## Loop over penalty functions and magnitudes


```r
policies <- sfSapply(penaltyfns, function(penalty) {
    policies <- sapply(c2, function(c2) {
        policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
            delta, reward, penalty = penalty(c2))
        i <- which(x_grid > K)[1]
        max(policycost$penalty_free_V[i, ])
    })
})
```


Note that `optim_policy` has been updated to return the equilibrium value of profits from fish harvests before the adjustment costs have been paid, `penalty_free_V`.  This containst the values for all possible states, we simply evaluate it at the carrying capacity (which is our initial condition.)  The index in `x_grid` that corresponds to the carrying capacity (initial condition) `i` indicates this.  



Quadratic costs on fishing effort have to be done separately,


```r
quad <- sapply(c2, function(c2) {
    effort_penalty = function(x, h) 0.1 * c2 * h/x
    policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
        delta, reward, penalty = fixed(0), effort_penalty)
    i <- which(x_grid > K)[1]
    max(policycost$penalty_free_V[i, ])  # chooses the most sensible harvest in t=1
})
dat <- cbind(policies, quad)
```


Tidy up the data and plot the net present value (before the penalty has been paid) relative to that achieved when managed without a penalty.  





























