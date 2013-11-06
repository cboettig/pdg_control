

```r
rm(list = ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```







```r
start <- Sys.time()
```



```r
compute_error_table <- function(reduction = 0.25, sigma = 0.2){
reduction <- reduction

price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)
c2 <- seq(0, 50, length.out=100)

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
ignore_fraction <- ignore_when_present$V1/optimal_cost$V1
assume_fraction <- assume_when_absent$V1/optimal_free$V1
assume_when_absent <- cbind(assume_when_absent, assume_fraction = assume_fraction)
ignore_when_present <- cbind(ignore_when_present, ignore_fraction = ignore_fraction)

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
    compute_error_table(r = r, sigma = s)    
  })
})
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
##    penalty_fn ignore_cost ignore_fraction assume_cost assume_fraction
## 1          L1     14536.1          0.9978       16858          1.0024
## 2          L2     11020.1          0.7584       15538          0.9240
## 3       fixed     13585.0          0.9207       17583          1.0455
## 4          L1      9273.7          1.0907       17562          1.0443
## 5          L2     11020.1          0.7584       15538          0.9240
## 6       fixed      8332.4          0.8351       17402          1.0348
## 7          L1      -641.8          0.2085       15451          0.9188
## 8          L2   -272586.1        -24.8951       11567          0.6878
## 9       fixed     -1213.0         -0.4593       15520          0.9229
## 10         L1     18282.0          1.0069       20974          0.9852
## 11         L2     12462.5          0.7202       19393          0.9110
## 12      fixed     16691.9          1.0269       20144          0.9463
## 13         L1     12789.0          0.9597       20663          0.9707
## 14         L2     12462.5          0.7202       19393          0.9110
## 15      fixed     12399.0          0.8846       21289          1.0001
## 16         L1      9099.5          0.7349       19999          0.9395
## 17         L2     -6337.8         -0.3916       17765          0.8345
## 18      fixed      7399.0          0.6234       21970          1.0320
## 19         L1     30227.4          1.0745       31377          0.9374
## 20         L2     33472.2          1.0000       33472          1.0000
## 21      fixed     28926.7          1.0890       30501          0.9112
## 22         L1     24157.9          1.1038       31768          0.9491
## 23         L2     19091.3          0.9305       29225          0.8731
## 24      fixed     23573.2          1.0250       30070          0.8984
## 25         L1     21086.1          1.1181       32036          0.9571
## 26         L2     19091.3          0.9305       29225          0.8731
## 27      fixed     17209.6          0.8533       32946          0.9843
##    sigma_g reduction L2 L1
## 1     0.05       0.1  1  1
## 2     0.05       0.1  1  1
## 3     0.05       0.1  1  1
## 4     0.05       0.2  2  1
## 5     0.05       0.2  2  1
## 6     0.05       0.2  2  1
## 7     0.05       0.3  3  1
## 8     0.05       0.3  3  1
## 9     0.05       0.3  3  1
## 10    0.20       0.1  1  2
## 11    0.20       0.1  1  2
## 12    0.20       0.1  1  2
## 13    0.20       0.2  2  2
## 14    0.20       0.2  2  2
## 15    0.20       0.2  2  2
## 16    0.20       0.3  3  2
## 17    0.20       0.3  3  2
## 18    0.20       0.3  3  2
## 19    0.50       0.1  1  3
## 20    0.50       0.1  1  3
## 21    0.50       0.1  1  3
## 22    0.50       0.2  2  3
## 23    0.50       0.2  2  3
## 24    0.50       0.2  2  3
## 25    0.50       0.3  3  3
## 26    0.50       0.3  3  3
## 27    0.50       0.3  3  3
```



```r
#df <- df[1:5]
library(xtable)
print(xtable(df), type="html")
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Tue Nov  5 08:20:47 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> penalty_fn </TH> <TH> ignore_cost </TH> <TH> ignore_fraction </TH> <TH> assume_cost </TH> <TH> assume_fraction </TH> <TH> sigma_g </TH> <TH> reduction </TH> <TH> L2 </TH> <TH> L1 </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> L1 </TD> <TD align="right"> 14536.09 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 16857.61 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> L2 </TD> <TD align="right"> 11020.10 </TD> <TD align="right"> 0.76 </TD> <TD align="right"> 15538.44 </TD> <TD align="right"> 0.92 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> fixed </TD> <TD align="right"> 13584.97 </TD> <TD align="right"> 0.92 </TD> <TD align="right"> 17582.78 </TD> <TD align="right"> 1.05 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> L1 </TD> <TD align="right"> 9273.66 </TD> <TD align="right"> 1.09 </TD> <TD align="right"> 17561.81 </TD> <TD align="right"> 1.04 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> L2 </TD> <TD align="right"> 11020.10 </TD> <TD align="right"> 0.76 </TD> <TD align="right"> 15538.44 </TD> <TD align="right"> 0.92 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> fixed </TD> <TD align="right"> 8332.45 </TD> <TD align="right"> 0.84 </TD> <TD align="right"> 17401.82 </TD> <TD align="right"> 1.03 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> L1 </TD> <TD align="right"> -641.78 </TD> <TD align="right"> 0.21 </TD> <TD align="right"> 15451.25 </TD> <TD align="right"> 0.92 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> L2 </TD> <TD align="right"> -272586.11 </TD> <TD align="right"> -24.90 </TD> <TD align="right"> 11567.12 </TD> <TD align="right"> 0.69 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> fixed </TD> <TD align="right"> -1213.01 </TD> <TD align="right"> -0.46 </TD> <TD align="right"> 15519.95 </TD> <TD align="right"> 0.92 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> L1 </TD> <TD align="right"> 18281.97 </TD> <TD align="right"> 1.01 </TD> <TD align="right"> 20973.56 </TD> <TD align="right"> 0.99 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> L2 </TD> <TD align="right"> 12462.51 </TD> <TD align="right"> 0.72 </TD> <TD align="right"> 19392.85 </TD> <TD align="right"> 0.91 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> fixed </TD> <TD align="right"> 16691.94 </TD> <TD align="right"> 1.03 </TD> <TD align="right"> 20143.80 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> L1 </TD> <TD align="right"> 12788.95 </TD> <TD align="right"> 0.96 </TD> <TD align="right"> 20663.21 </TD> <TD align="right"> 0.97 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> L2 </TD> <TD align="right"> 12462.51 </TD> <TD align="right"> 0.72 </TD> <TD align="right"> 19392.85 </TD> <TD align="right"> 0.91 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> fixed </TD> <TD align="right"> 12399.01 </TD> <TD align="right"> 0.88 </TD> <TD align="right"> 21289.21 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> L1 </TD> <TD align="right"> 9099.48 </TD> <TD align="right"> 0.73 </TD> <TD align="right"> 19999.07 </TD> <TD align="right"> 0.94 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> L2 </TD> <TD align="right"> -6337.76 </TD> <TD align="right"> -0.39 </TD> <TD align="right"> 17765.33 </TD> <TD align="right"> 0.83 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> fixed </TD> <TD align="right"> 7399.01 </TD> <TD align="right"> 0.62 </TD> <TD align="right"> 21969.84 </TD> <TD align="right"> 1.03 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> L1 </TD> <TD align="right"> 30227.43 </TD> <TD align="right"> 1.07 </TD> <TD align="right"> 31377.07 </TD> <TD align="right"> 0.94 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> L2 </TD> <TD align="right"> 33472.18 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 33472.18 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> fixed </TD> <TD align="right"> 28926.73 </TD> <TD align="right"> 1.09 </TD> <TD align="right"> 30500.85 </TD> <TD align="right"> 0.91 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> L1 </TD> <TD align="right"> 24157.93 </TD> <TD align="right"> 1.10 </TD> <TD align="right"> 31767.92 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> L2 </TD> <TD align="right"> 19091.31 </TD> <TD align="right"> 0.93 </TD> <TD align="right"> 29225.46 </TD> <TD align="right"> 0.87 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> fixed </TD> <TD align="right"> 23573.19 </TD> <TD align="right"> 1.02 </TD> <TD align="right"> 30069.80 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> L1 </TD> <TD align="right"> 21086.12 </TD> <TD align="right"> 1.12 </TD> <TD align="right"> 32035.93 </TD> <TD align="right"> 0.96 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> L2 </TD> <TD align="right"> 19091.31 </TD> <TD align="right"> 0.93 </TD> <TD align="right"> 29225.46 </TD> <TD align="right"> 0.87 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> fixed </TD> <TD align="right"> 17209.56 </TD> <TD align="right"> 0.85 </TD> <TD align="right"> 32946.24 </TD> <TD align="right"> 0.98 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
   </TABLE>

```r
writeLines(print(xtable(df)), "output.md")
```

% latex table generated in R 3.0.2 by xtable 1.7-1 package
% Tue Nov  5 08:20:47 2013
\begin{table}[ht]
\centering
\begin{tabular}{rlrrrrrrrr}
  \hline
 & penalty\_fn & ignore\_cost & ignore\_fraction & assume\_cost & assume\_fraction & sigma\_g & reduction & L2 & L1 \\ 
  \hline
1 & L1 & 14536.09 & 1.00 & 16857.61 & 1.00 & 0.05 & 0.10 &   1 &   1 \\ 
  2 & L2 & 11020.10 & 0.76 & 15538.44 & 0.92 & 0.05 & 0.10 &   1 &   1 \\ 
  3 & fixed & 13584.97 & 0.92 & 17582.78 & 1.05 & 0.05 & 0.10 &   1 &   1 \\ 
  4 & L1 & 9273.66 & 1.09 & 17561.81 & 1.04 & 0.05 & 0.20 &   2 &   1 \\ 
  5 & L2 & 11020.10 & 0.76 & 15538.44 & 0.92 & 0.05 & 0.20 &   2 &   1 \\ 
  6 & fixed & 8332.45 & 0.84 & 17401.82 & 1.03 & 0.05 & 0.20 &   2 &   1 \\ 
  7 & L1 & -641.78 & 0.21 & 15451.25 & 0.92 & 0.05 & 0.30 &   3 &   1 \\ 
  8 & L2 & -272586.11 & -24.90 & 11567.12 & 0.69 & 0.05 & 0.30 &   3 &   1 \\ 
  9 & fixed & -1213.01 & -0.46 & 15519.95 & 0.92 & 0.05 & 0.30 &   3 &   1 \\ 
  10 & L1 & 18281.97 & 1.01 & 20973.56 & 0.99 & 0.20 & 0.10 &   1 &   2 \\ 
  11 & L2 & 12462.51 & 0.72 & 19392.85 & 0.91 & 0.20 & 0.10 &   1 &   2 \\ 
  12 & fixed & 16691.94 & 1.03 & 20143.80 & 0.95 & 0.20 & 0.10 &   1 &   2 \\ 
  13 & L1 & 12788.95 & 0.96 & 20663.21 & 0.97 & 0.20 & 0.20 &   2 &   2 \\ 
  14 & L2 & 12462.51 & 0.72 & 19392.85 & 0.91 & 0.20 & 0.20 &   2 &   2 \\ 
  15 & fixed & 12399.01 & 0.88 & 21289.21 & 1.00 & 0.20 & 0.20 &   2 &   2 \\ 
  16 & L1 & 9099.48 & 0.73 & 19999.07 & 0.94 & 0.20 & 0.30 &   3 &   2 \\ 
  17 & L2 & -6337.76 & -0.39 & 17765.33 & 0.83 & 0.20 & 0.30 &   3 &   2 \\ 
  18 & fixed & 7399.01 & 0.62 & 21969.84 & 1.03 & 0.20 & 0.30 &   3 &   2 \\ 
  19 & L1 & 30227.43 & 1.07 & 31377.07 & 0.94 & 0.50 & 0.10 &   1 &   3 \\ 
  20 & L2 & 33472.18 & 1.00 & 33472.18 & 1.00 & 0.50 & 0.10 &   1 &   3 \\ 
  21 & fixed & 28926.73 & 1.09 & 30500.85 & 0.91 & 0.50 & 0.10 &   1 &   3 \\ 
  22 & L1 & 24157.93 & 1.10 & 31767.92 & 0.95 & 0.50 & 0.20 &   2 &   3 \\ 
  23 & L2 & 19091.31 & 0.93 & 29225.46 & 0.87 & 0.50 & 0.20 &   2 &   3 \\ 
  24 & fixed & 23573.19 & 1.02 & 30069.80 & 0.90 & 0.50 & 0.20 &   2 &   3 \\ 
  25 & L1 & 21086.12 & 1.12 & 32035.93 & 0.96 & 0.50 & 0.30 &   3 &   3 \\ 
  26 & L2 & 19091.31 & 0.93 & 29225.46 & 0.87 & 0.50 & 0.30 &   3 &   3 \\ 
  27 & fixed & 17209.56 & 0.85 & 32946.24 & 0.98 & 0.50 & 0.30 &   3 &   3 \\ 
   \hline
\end{tabular}
\end{table}

```r
system("pandoc output.md -o output.pdf")
write.csv(df, "output.csv")
```



```r
Sys.time() - start 
```

```
## Time difference of 1.01 hours
```

