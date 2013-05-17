

```r
rm(list = ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```






Define the parameters of our cost function. We will 
consider a fixed price and fixed $c_0$ cost per unit effort.  
An additional "policy cost" is introduced through the $c_2$
coefficient, which can take a range of values. We will loop
over this range for each of the functional forms of the policy
cost, in order to choose coefficients $c_2$ in each case
that are comparable.  



```r
price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)
c2 <- exp(seq(0, log(21), length.out = 20)) - 1
```











This block defines the various parameters and nuisance parameters we
need to specify our optimal control problem.



```r
seed <- 123  # Random seed (replicable results)
delta <- 0.05  # economic discounting rate
OptTime <- 20  # stopping time
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
x_grid <- seq(0.01, 1.2 * K, length = gridsize)
h_grid <- seq(0.01, 0.8 * K, length = gridsize)
```



Given these parameters, we can determine the optimal solution under the
classic assumption of "adjustment-free" costs, where there is no penalty
for adjusting the value of the control (harvest) at each decision step
(year).



```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward = reward)
```



Now we introduce the functional forms for each of the adjustment costs,



```r
L1 <- function(c2) function(h, h_prev) c2 * abs(h - h_prev)
free_increase <- function(c2) function(h, h_prev) c2 * abs(min(h - h_prev, 0))  # increasing harvest is free
free_decrease <- function(c2) function(h, h_prev) c2 * max(h - h_prev, 0)  # decreasing harvest is free
fixed <- function(c2) function(h, h_prev) c2 * as.numeric(!(h == h_prev))
L2 <- function(c2) function(h, h_prev) c2 * (h - h_prev)^2
none <- function(h, h_prev) 0
penaltyfns <- list(L2 = L2, L1 = L1, free_decrease = free_decrease, fixed = fixed, 
    free_increase = free_increase)
```



## Apples to Apples levels

Before we can start comparing solutions induced under each of these
costs, we need to scale them appropriately (e.g. the coefficient $c_2$
has different units and a different magnitude of effect when multiplied
by quadratic differences $(h_{t} - h_{t-1})^2$ then when multiplied by
linear differences $h_{t} - h_{t-1}$.


To do so we loop over a grid of $c_2$ values and determine the net
present value under each functional form of adjustment costs for each
$c_2$ in our grid.



```r
require(snowfall)
sfInit(cpu = 16, parallel = T)
```

```
## R Version:  R version 2.15.2 (2012-10-26)
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 16
## CPUs.
```

```r
sfLibrary(pdgControl)
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Warning: 'keep.source' is deprecated and will be ignored
```

```r
sfExportAll()
```



### Loop over penalty functions and magnitudes


```r
policies <- lapply(penaltyfns, function(penalty) {
    sfLapply(c2, function(c2) {
        policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
            delta, reward, penalty = penalty(c2))
    })
})
```


Note that `optim_policy` has been updated to return the equilibrium value
of profits from fish harvests before the adjustment costs have been paid,
`penalty_free_V`, as well as the total value `V`. Initially we based the 
comparison on matching `penalty_free_V`, now the economists recommend we
simply use the total net present value resulting, `V`.  



The value matrices returned in either case contain the values for all possible states,
we simply evaluate it at the carrying capacity (which is our initial
condition.)  The index in `x_grid` that corresponds to the carrying
capacity (initial condition) `i` indicates this.


We also do the the classical quadratic costs on fishing effort separately,


```r
quad <- sfLapply(c2, function(c2) {
    effort_penalty = function(x, h) (c2 * h/x)/price
    policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
        delta, reward, penalty = none, effort_penalty)
})
policies <- c(policies, quad = list(quad))
```



Extract the policy cost 


```r
i <- which(x_grid > K)[1]
fees <- lapply(policies, function(penalty) sapply(penalty, function(c2_run) max(c2_run$V[i, 
    ])  # Would be penalty_free_V originally
))
```



Tidy up the data and plot the net present value (before the penalty has
been paid) relative to that achieved when managed without a penalty.


```r
npv0 <- max(fees$L1)  # all have same max, at c2=0
npv0
```

```
## [1] 183
```

```r
fees <- data.frame(c2 = c2, fees)
fees <- melt(fees, id = "c2")
ggplot(fees, aes(c2, value, col = variable)) + geom_point() + geom_line()
```

![plot of chunk npv-plot](http://farm8.staticflickr.com/7318/8748417848_fecba4020e_o.png) 


Find the value of `c2` that brings each penalty closest to 75% of the cost-free adjustment value:


```r
ggplot(fees, aes(c2, (npv0 - value)/npv0, col = variable)) + geom_point() + 
    geom_line()
```

![plot of chunk apples_plot](http://farm8.staticflickr.com/7307/8748417934_26cc1eb952_o.png) 









```r
closest <- function(x, v) {
    which.min(abs(v - x))
}
dt_npv <- data.table(fees)
index <- dt_npv[, closest(0.25, (npv0 - value)/npv0), by = variable]
apples_index <- index$V1
names(apples_index) = index$variable
apples <- c2[index$V1]
names(apples) = index$variable
apples
```

```
##            L2            L1 free_decrease         fixed free_increase 
##         1.228         2.070         0.000        10.063         0.000 
##          quad 
##         2.603
```


Solve the policy cost for the specified penalty function


```r
L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D
free_increase_policy <- policies$free_increase[[apples_index["free_increase"]]]$D
free_decrease_policy <- policies$free_decrease[[apples_index["free_decrease"]]]$D
quad_policy <- policies$quad[[apples_index["quad"]]]$D
```



```r
quad_profit <- profit_harvest(price = price, c0 = c0, c1 = apples["quad"])
sims <- list(L1 = simulate_optim(f, pars, x_grid, h_grid, x0, L1_policy, z_g, 
    z_m, z_i, opt$D, profit = profit, penalty = L1(apples["L1"]), seed = seed), 
    L2 = simulate_optim(f, pars, x_grid, h_grid, x0, L2_policy, z_g, z_m, z_i, 
        opt$D, profit = profit, penalty = L2(apples["L2"]), seed = seed), fixed = simulate_optim(f, 
        pars, x_grid, h_grid, x0, fixed_policy, z_g, z_m, z_i, opt$D, profit = profit, 
        penalty = fixed(apples["fixed"]), seed = seed), increase = simulate_optim(f, 
        pars, x_grid, h_grid, x0, free_increase_policy, z_g, z_m, z_i, opt$D, 
        profit = profit, penalty = free_increase(apples["increase"]), seed = seed), 
    decrease = simulate_optim(f, pars, x_grid, h_grid, x0, free_decrease_policy, 
        z_g, z_m, z_i, opt$D, profit = profit, penalty = free_decrease(apples["decrease"]), 
        seed = seed), quad = simulate_optim(f, pars, x_grid, h_grid, x0, quad_policy, 
        z_g, z_m, z_i, opt$D, profit = quad_profit, penalty = none, seed = seed))
```




```r
# Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id = names(sims[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "penalty_fn")  # names are nice
```





```r
v <- dt[, var(harvest), by = "penalty_fn"]
var <- v$V1
names(var) <- v$penalty_fn
acor <- dt[, acf(harvest, plot = F)$acf[2], by = "penalty_fn"]$V1
names(acor) <- names(var)
out <- rbind(var = var, a = acor)
out
```

```
##          L1     L2   fixed increase decrease    quad
## var  3.5800 0.8466  8.5305   6.1348   6.1348  1.6863
## a   -0.1617 0.2210 -0.2585  -0.3741  -0.3741 -0.1726
```






# Plots 



```r
p1 <- ggplot(dt) + geom_line(aes(time, alternate), col = "grey20", lwd = 1) + 
    geom_line(aes(time, fishstock), col = rgb(0, 0, 1, 0.8)) + facet_wrap(~penalty_fn) + 
    labs(x = "time", y = "stock size", title = "Stock Dynamics")
p1
```

![plot of chunk p1](http://farm8.staticflickr.com/7318/8747296307_01e8270259_o.png) 



```r
p2 <- ggplot(dt) + geom_line(aes(time, harvest_alt), col = "grey20", lwd = 1) + 
    geom_line(aes(time, harvest), col = rgb(0, 0, 1, 0.8)) + facet_wrap(~penalty_fn) + 
    labs(x = "time", y = "havest intensity (fish taken)", title = "Harvest Policy Dynamics")
p2
```

![plot of chunk Figure3](http://farm8.staticflickr.com/7322/8747296383_6f43b33cef_o.png) 


# Figure 4


```r
frac_lost <- seq(0, 1, length = 20)
```



```r
fig4 <- function(fraction_lost) {
    closest <- function(x, v) {
        which.min(abs(v - x))
    }
    dt_npv <- data.table(fees)
    index <- dt_npv[, closest(fraction_lost, (npv0 - value)/npv0), by = variable]
    apples_index <- index$V1
    names(apples_index) = index$variable
    apples <- c2[index$V1]
    names(apples) = index$variable
    
    L2_policy <- policies$L2[[apples_index["L2"]]]$D
    L1_policy <- policies$L1[[apples_index["L1"]]]$D
    fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D
    free_increase_policy <- policies$free_increase[[apples_index["free_increase"]]]$D
    free_decrease_policy <- policies$free_decrease[[apples_index["free_decrease"]]]$D
    quad_policy <- policies$quad[[apples_index["quad"]]]$D
    
    quad_profit <- profit_harvest(price = price, c0 = c0, c1 = apples["quad"])
    sims <- lapply(1:50, function(reps) list(L1 = simulate_optim(f, pars, x_grid, 
        h_grid, x0, L1_policy, z_g, z_m, z_i, opt$D, profit = profit, penalty = L1(apples["L1"])), 
        L2 = simulate_optim(f, pars, x_grid, h_grid, x0, L2_policy, z_g, z_m, 
            z_i, opt$D, profit = profit, penalty = L2(apples["L2"])), fixed = simulate_optim(f, 
            pars, x_grid, h_grid, x0, fixed_policy, z_g, z_m, z_i, opt$D, profit = profit, 
            penalty = fixed(apples["fixed"])), increase = simulate_optim(f, 
            pars, x_grid, h_grid, x0, free_increase_policy, z_g, z_m, z_i, opt$D, 
            profit = profit, penalty = free_increase(apples["increase"])), decrease = simulate_optim(f, 
            pars, x_grid, h_grid, x0, free_decrease_policy, z_g, z_m, z_i, opt$D, 
            profit = profit, penalty = free_decrease(apples["decrease"])), quad = simulate_optim(f, 
            pars, x_grid, h_grid, x0, quad_policy, z_g, z_m, z_i, opt$D, profit = quad_profit, 
            penalty = none)))
    
    # Make data tidy (melt), fast (data.tables), and nicely labeled.
    dat <- melt(sims, id = names(sims[[1]][[1]]))
    dt <- data.table(dat)
    setnames(dt, "L1", "reps")  # names are nice
    setnames(dt, "L2", "penalty_fn")  # names are nice
    
    dt
}
```




```r
# Helper functions to extract the summary stats in different variables
stats_harvest <- function(dt) {
    v <- dt[, var(harvest), by = c("penalty_fn", "reps")]
    acor <- dt[, acf(harvest, plot = F)$acf[2], by = c("penalty_fn", "reps")]
    df <- cbind(v, acor$V1)
    setnames(df, c("V1", "V2"), c("var", "acor"))  # names are nice
    df
}

stats_fishstock <- function(dt) {
    v <- dt[, var(fishstock), by = c("penalty_fn", "reps")]
    acor <- dt[, acf(fishstock, plot = F)$acf[2], by = c("penalty_fn", "reps")]
    df <- cbind(v, acor$V1)
    setnames(df, c("V1", "V2"), c("var", "acor"))  # names are nice
    df
}
```




```r
#' a simple function for reorganizing the data over the different 'faction-lost' levels
get_trends <- function(tmp2) {
    out <- melt(tmp2, id = c("reps", "var", "acor"))
    colnames(out) = c("reps", "var", "acor", "nothing", "penalty", "index")
    out <- cbind(out[c(1, 2, 3, 5)], fraction = frac_lost[out$index])
    out <- data.table(out)
    Ev = out[, mean(var), by = c("penalty", "fraction")]
    SDv = out[, sd(var), by = c("penalty", "fraction")]
    Ea = out[, mean(acor), by = c("penalty", "fraction")]
    SDa = out[, sd(acor), by = c("penalty", "fraction")]
    harvest_trends <- data.table(penalty = Ev$penalty, fraction = Ev$fraction, 
        Ev = Ev$V1, Ea = Ea$V1, SDv = SDv$V1, SDa = SDa$V1)
}
```






```r
sims_at_each_apple <- lapply(frac_lost, fig4)
harvest_stats <- lapply(sims_at_each_apple, stats_harvest)
harvest_trends <- get_trends(harvest_stats)
```



```r
Fig4a <- ggplot(harvest_trends, aes(fraction, Ev, ymin = Ev - SDv, ymax = Ev + 
    SDv, col = penalty)) + geom_ribbon(aes(fill = penalty, col = NA), lwd = 0, 
    alpha = 0.05) + geom_line() + xlab("Fraction of NPV lost to costs")

Fig4b <- ggplot(harvest_trends, aes(fraction, Ea, ymin = Ea - SDa, ymax = Ea + 
    SDa, col = penalty)) + geom_ribbon(aes(fill = penalty, col = NA), lwd = 0, 
    alpha = 0.05) + geom_line() + xlab("Fraction of NPV lost to costs")
Fig4a
```

![plot of chunk Figure4](http://farm9.staticflickr.com/8134/8747296495_bbdb3bd9b9_o.png) 

```r
Fig4b
```

![plot of chunk Figure4](http://farm8.staticflickr.com/7304/8747296589_857eb8a5a5_o.png) 




```r
fishstock_stats <- lapply(sims_at_each_apple, stats_fishstock)
fishstock_trends <- get_trends(fishstock_stats)

FigS4a <- ggplot(fishstock_trends, aes(fraction, Ev, ymin = Ev - SDv, ymax = Ev + 
    SDv, col = penalty)) + geom_ribbon(aes(fill = penalty, col = NA), lwd = 0, 
    alpha = 0.05) + geom_line() + xlab("Fraction of NPV lost to costs")

FigS4b <- ggplot(fishstock_trends, aes(fraction, Ea, ymin = Ea - SDa, ymax = Ea + 
    SDa, col = penalty)) + geom_ribbon(aes(fill = penalty, col = NA), lwd = 0, 
    alpha = 0.05) + geom_line() + xlab("Fraction of NPV lost to costs")
FigS4a
```

![plot of chunk Figure4S](http://farm9.staticflickr.com/8544/8747296671_08283b6117_o.png) 

```r
FigS4b
```

![plot of chunk Figure4S](http://farm9.staticflickr.com/8552/8747296745_1c2f43bc42_o.png) 












