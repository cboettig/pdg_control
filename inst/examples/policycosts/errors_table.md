

```r
rm(list = ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```









```r
compute_error_table <- function(reduction = 0.25, sigma = 0.2){
reduction <- reduction

price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)
c2 <- exp(seq(0, log(41), length.out = 40))-1
c2 <- seq(0, 40, length.out=100)

seed <- 123                 # Random seed (replicable results)
delta <- 0.05               # economic discounting rate
OptTime <- 20               # stopping time
gridsize <- 50              # grid size for fish stock and harvest rate (discretized population)
sigma_g <- sigma              # Noise in population growth
reward <- 0                 # bonus for satisfying the boundary condition
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() 1         # No measurement noise, 
z_i <- function() 1         # No implemenation noise
f <- BevHolt                # Select the state equation
pars <- c(1.5, 0.05)        # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference 
xT <- 0                     # boundary conditions
x0 <- K
x_grid <- seq(0.01, 1.2 * K, length = gridsize)  
h_grid <- seq(0.01, 0.8 * K, length = gridsize)  

SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)

L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
free_increase <- function(c2) function(h, h_prev)  c2 * abs(pmin(h - h_prev, 0)) # increasing harvest is free
free_decrease <- function(c2) function(h, h_prev)  c2 * pmax(h - h_prev, 0) # decreasing harvest is free
fixed <-  function(c2) function(h, h_prev) c2 * as.numeric( !(h == h_prev) )
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2
none <- function(h, h_prev)  0
penaltyfns <- list(L2=L2, L1=L1, free_decrease=free_decrease, fixed=fixed, free_increase=free_increase)

require(snowfall)
sfInit(cpu=8, parallel=T)
sfLibrary(pdgControl)
sfExportAll()

policies <- lapply(penaltyfns, function(penalty){
  sfLapply(c2, function(c2){
      policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                   profit, delta, reward, penalty = penalty(c2))
      }
  )
})

quad <- 
  sfLapply(c2, function(c2){
  effort_penalty = function(x,h) (c2 * h / x) / price
  policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                        profit, delta, reward, penalty = none, 
                        effort_penalty)
})
policies <- c(policies, quad=list(quad))

i <- which(x_grid > K)[1]
fees <- 
lapply(policies, function(penalty) 
  sapply(penalty, function(c2_run)
    max(c2_run$V[i,]) # Would be penalty_free_V originally   ## this isn't correct for asym cases
  )
)

npv0 <- max(fees$L1) # all have same max, at c2=0 
npv0
fees <- data.frame(c2=c2,fees)
fees <- melt(fees, id="c2")
#ggplot(fees, aes(c2, value, col=variable)) + geom_point() + geom_line()

closest <- function(x, v){
  which.min(abs(v-x))
}
dt_npv <- data.table(fees)
index <- dt_npv[,closest(reduction, (npv0-value)/npv0), by=variable]
apples_index <- index$V1
names(apples_index) = index$variable
apples <- c2[index$V1]
names(apples) = index$variable
apples

setkey(dt_npv, variable)
values <- apply(index, 1, function(x) dt_npv[x[1], ][as.integer(x[2]), ]$value)

percent.error <- (values - ((1-reduction)*npv0)) / ((1-reduction)*npv0)* 100

print_npv <- data.frame(model=index$variable, "value realized"=values, 
                        "percent of npv0" = 100*values/npv0, "percent error"=percent.error)


L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D
free_increase_policy <- policies$free_increase[[apples_index["free_increase"]]]$D
free_decrease_policy <- policies$free_decrease[[apples_index["free_decrease"]]]$D
quad_policy <- policies$quad[[apples_index["quad"]]]$D

quad_profit <- profit_harvest(price = price, c0 = c0, c1 = apples["quad"]) 

reps <- 1:100
names(reps) = paste("rep", 1:10, sep="_")
sims <- list(
  L1 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L1_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L1(apples["L1"]), seed=seed)), 
  L2 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L2_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L2(apples["L2"]), seed=seed)),
  fixed = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                         fixed_policy, z_g, z_m, z_i, 
                         opt$D, profit=profit, penalty=fixed(apples["fixed"]), seed=seed))
  
)

#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]][[1]]))  
dt <- data.table(dat)
setnames(dt, "L2", "replicate") # names are nice
setnames(dt, "L1", "penalty_fn") # names are nice

sims_free <- list(
  L1 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L1_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=none, seed=seed)), 
  L2 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L2_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=none, seed=seed)),
  fixed = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                         fixed_policy, z_g, z_m, z_i, 
                         opt$D, profit=profit, penalty=none, seed=seed))
)

#Make data tidy (melt), fast (data.tables), and nicely labeled.
free <- melt(sims_free, id=names(sims_free[[1]][[1]]))  
dt_free <- data.table(free)
setnames(dt_free, "L2", "replicate") # names are nice
setnames(dt_free, "L1", "penalty_fn") # names are nice

# Profit when accounting for penalty when present
optimal_cost <- dt[, sum(profit_fishing - policy_cost), by=penalty_fn ] 
# Profit when ignoring penalty when present
ignore_when_present <- dt[, sum(profit_fishing_alt - policy_cost_alt), by=penalty_fn] 
# Profit when assuming penalty when it is absent
assume_when_absent <- dt_free[, sum(profit_fishing - policy_cost), by=penalty_fn]
# Profit when ignoring penalty when it is absent
optimal_free <- dt_free[, sum(profit_fishing_alt - policy_cost_alt), by=penalty_fn]

# Normalize by the optimal 
#ignore_when_present$V1 <- ignore_when_present$V1/optimal_cost$V1
#assume_when_absent$V1 <- assume_when_absent$V1/optimal_free$V1

# Name and merge columns
setnames(ignore_when_present, "V1", "ignore_cost")
setnames(assume_when_absent, "V1", "assume_cost")
error_costs = merge(ignore_when_present, assume_when_absent, "penalty_fn")

print_npv # theoretically acheivable profits under these costs
optimal_cost$V1/optimal_free$V1 # actually realized 

error_costs <- cbind(error_costs, sigma_g = sigma_g, reduction = reduction)

error_costs 
}
```



```r
sigmas <- c(0.05, 0.2, 0.5)
reductions <- c(0.1, 0.2, 0.3)
out <- lapply(sigmas, function(s){
  lapply(reductions, function(r){
    compute_error_table(r = 0.25, sigma = s)    
  })
})
```

```
## Loading required package: snowfall
```

```
## Loading required package: snow
```

```
## R Version:  R version 3.0.2 (2013-09-25)
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```
## Explicit sfStop() is missing: stop now.
```

```
## Stopping cluster
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 8
## CPUs.
```

```
## Library pdgControl loaded.
```

```
## Library pdgControl loaded in cluster.
```

```r
       
who <- names(out[[1]][[1]])
df <- melt(out, id=who)              
df              
```

```
##    penalty_fn ignore_cost assume_cost sigma_g reduction L2 L1
## 1          L1        4800       17576    0.05      0.25  1  1
## 2          L2       -5883       14153    0.05      0.25  1  1
## 3       fixed        5787       15520    0.05      0.25  1  1
## 4          L1        4800       17576    0.05      0.25  2  1
## 5          L2       -5883       14153    0.05      0.25  2  1
## 6       fixed        5787       15520    0.05      0.25  2  1
## 7          L1        4800       17576    0.05      0.25  3  1
## 8          L2       -5883       14153    0.05      0.25  3  1
## 9       fixed        5787       15520    0.05      0.25  3  1
## 10         L1       10450       18858    0.20      0.25  1  2
## 11         L2        3565       18214    0.20      0.25  1  2
## 12      fixed       10177       21010    0.20      0.25  1  2
## 13         L1       10450       18858    0.20      0.25  2  2
## 14         L2        3565       18214    0.20      0.25  2  2
## 15      fixed       10177       21010    0.20      0.25  2  2
## 16         L1       10450       18858    0.20      0.25  3  2
## 17         L2        3565       18214    0.20      0.25  3  2
## 18      fixed       10177       21010    0.20      0.25  3  2
## 19         L1       23537       31768    0.50      0.25  1  3
## 20         L2       21416       29918    0.50      0.25  1  3
## 21      fixed       20462       33782    0.50      0.25  1  3
## 22         L1       23537       31768    0.50      0.25  2  3
## 23         L2       21416       29918    0.50      0.25  2  3
## 24      fixed       20462       33782    0.50      0.25  2  3
## 25         L1       23537       31768    0.50      0.25  3  3
## 26         L2       21416       29918    0.50      0.25  3  3
## 27      fixed       20462       33782    0.50      0.25  3  3
```

