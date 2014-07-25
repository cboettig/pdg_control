

```r
rm(list=ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```





```r
start <- Sys.time()
```


```r
compute_error_table <- function(r, sigma){
reduction_range = c(.1, .2, .3)
reduction <- reduction_range[r]

price = 10
c0 = 0
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)
c2 <- seq(0, 50, length.out=100)


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
fixed <-  function(c2) function(h, h_prev) c2 * as.numeric( !(h == h_prev) )
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2
none <- function(h, h_prev)  0


reduction_table <- read.csv('reduction_table.csv', row.names=1)
apples <-   reduction_table[[r]]
names(apples) <- rownames(reduction_table)

L2_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                   profit, delta, reward, penalty = L2(apples[["L2"]]))
L1_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                   profit, delta, reward, penalty = L1(apples[["L1"]]))
fixed_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                   profit, delta, reward, penalty = fixed(apples[["fixed"]]))


reps <- 1:100
names(reps) = paste("rep", 1:100, sep="_")
seeds <- 1:100
sims <- list(
  L1 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L1_policy$D, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L1(apples[["L1"]]), seed=seeds[x])), 
  L2 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L2_policy$D, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L2(apples[["L2"]]), seed=seeds[x])),
  fixed = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                         fixed_policy$D, z_g, z_m, z_i, 
                         opt$D, profit=profit, penalty=fixed(apples[["fixed"]]), seed=seeds[x]))
  
)

#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]][[1]]))  
dt <- data.table(dat)
setnames(dt, "L2", "replicate") # names are nice
setnames(dt, "L1", "penalty_fn") # names are nice

dt$profit <- dt$profit_fishing/(1+delta)^dt$time
dt$profit_fishing_alt <- dt$profit_fishing_alt/(1+delta)^dt$time
dt$policy_cost <- dt$policy_cost/(1+delta)^dt$time
dt$policy_cost_alt <- dt$policy_cost_alt/(1+delta)^dt$time

# Profit when accounting for penalty when present
optimal_cost <- dt[, sum(profit_fishing - policy_cost), by=penalty_fn ] 
# Profit when ignoring penalty when present
ignore_when_present <- dt[, sum(profit_fishing_alt - policy_cost_alt), by=penalty_fn] 
# Profit when assuming penalty when it is absent
assume_when_absent <- dt[, sum(profit_fishing), by=penalty_fn]
# Profit when ignoring penalty when it is absent
optimal_free <- dt[, sum(profit_fishing_alt), by=penalty_fn]

# Normalize by the optimal 
ignore_fraction <- ignore_when_present$V1/optimal_free$V1 # common normalization
assume_fraction <- assume_when_absent$V1/optimal_free$V1
assume_when_absent <- cbind(assume_when_absent, assume_fraction = assume_fraction, normalize_optimal_free=optimal_free$V1, normalize_optimal_cost = optimal_cost$V1)
ignore_when_present <- cbind(ignore_when_present, ignore_fraction = ignore_fraction)

# Name and merge columns
setnames(ignore_when_present, "V1", "ignore_cost")
setnames(assume_when_absent, "V1", "assume_cost")
error_costs = merge(ignore_when_present, assume_when_absent, "penalty_fn")


optimal_cost$V1/optimal_free$V1 # actually realized 

error_costs <- cbind(error_costs, sigma_g = sigma_g, reduction = reduction)

error_costs 
}
```


```r
sigmas <- c(0.01, 0.1, 0.2)
reductions <- 1:3
out <- lapply(sigmas, function(s){
  lapply(reductions, function(r){
    compute_error_table(r = r, sigma = s)    
  })
})
       
who <- names(out[[1]][[1]])
df <- melt(out, id=who)              
df
```

```
##    penalty_fn ignore_cost ignore_fraction assume_cost assume_fraction
## 1          L1       17188          0.9353       28052          1.5265
## 2          L2       18377          1.0000       28096          1.5289
## 3       fixed       15285          0.8317       27761          1.5107
## 4          L1       16593          0.9029       28052          1.5265
## 5          L2       16632          0.9050       27769          1.5111
## 6       fixed       11423          0.6216       27029          1.4708
## 7          L1       15411          0.8386       27970          1.5221
## 8          L2       15305          0.8328       27396          1.4908
## 9       fixed        4216          0.2294       10796          0.5875
## 10         L1       17159          0.9277       28442          1.5378
## 11         L2       18495          1.0000       28531          1.5426
## 12      fixed       14887          0.8049       27808          1.5035
## 13         L1       16503          0.8923       28413          1.5363
## 14         L2       16702          0.9030       28176          1.5234
## 15      fixed       10457          0.5654       27110          1.4658
## 16         L1       15192          0.8214       28340          1.5323
## 17         L2       15193          0.8215       27861          1.5064
## 18      fixed        4335          0.2344       10783          0.5830
## 19         L1       17659          0.9160       29755          1.5433
## 20         L2       19280          1.0000       29885          1.5501
## 21      fixed       15761          0.8175       29271          1.5182
## 22         L1       16872          0.8751       29675          1.5392
## 23         L2       17197          0.8920       29469          1.5285
## 24      fixed       11466          0.5947       28761          1.4918
## 25         L1       15312          0.7942       29516          1.5309
## 26         L2       15477          0.8028       29060          1.5073
## 27      fixed        6761          0.3507       12609          0.6540
##    normalize_optimal_free normalize_optimal_cost sigma_g reduction L2 L1
## 1                   18377                  26898    0.01       0.1  1  1
## 2                   18377                  28096    0.01       0.1  1  1
## 3                   18377                  25901    0.01       0.1  1  1
## 4                   18377                  26322    0.01       0.2  2  1
## 5                   18377                  26919    0.01       0.2  2  1
## 6                   18377                  22992    0.01       0.2  2  1
## 7                   18377                  25154    0.01       0.3  3  1
## 8                   18377                  26559    0.01       0.3  3  1
## 9                   18377                   9681    0.01       0.3  3  1
## 10                  18495                  27334    0.10       0.1  1  2
## 11                  18495                  28531    0.10       0.1  1  2
## 12                  18495                  25673    0.10       0.1  1  2
## 13                  18495                  26802    0.10       0.2  2  2
## 14                  18495                  27451    0.10       0.2  2  2
## 15                  18495                  23144    0.10       0.2  2  2
## 16                  18495                  25753    0.10       0.3  3  2
## 17                  18495                  26967    0.10       0.3  3  2
## 18                  18495                   9667    0.10       0.3  3  2
## 19                  19280                  28390    0.20       0.1  1  3
## 20                  19280                  29885    0.20       0.1  1  3
## 21                  19280                  27189    0.20       0.1  1  3
## 22                  19280                  27770    0.20       0.2  2  3
## 23                  19280                  28396    0.20       0.2  2  3
## 24                  19280                  24712    0.20       0.2  2  3
## 25                  19280                  26520    0.20       0.3  3  3
## 26                  19280                  27709    0.20       0.3  3  3
## 27                  19280                  10927    0.20       0.3  3  3
```


```r
#df <- df[1:5]
library(xtable)
print(xtable(df), type="html")
```

<!-- html table generated in R 3.1.0 by xtable 1.7-3 package -->
<!-- Tue Jun  3 15:05:26 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> penalty_fn </TH> <TH> ignore_cost </TH> <TH> ignore_fraction </TH> <TH> assume_cost </TH> <TH> assume_fraction </TH> <TH> normalize_optimal_free </TH> <TH> normalize_optimal_cost </TH> <TH> sigma_g </TH> <TH> reduction </TH> <TH> L2 </TH> <TH> L1 </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> L1 </TD> <TD align="right"> 17187.64 </TD> <TD align="right"> 0.94 </TD> <TD align="right"> 28052.13 </TD> <TD align="right"> 1.53 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 26898.38 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> L2 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 28095.87 </TD> <TD align="right"> 1.53 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 28095.87 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> fixed </TD> <TD align="right"> 15284.62 </TD> <TD align="right"> 0.83 </TD> <TD align="right"> 27760.96 </TD> <TD align="right"> 1.51 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 25900.76 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> L1 </TD> <TD align="right"> 16593.12 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 28052.13 </TD> <TD align="right"> 1.53 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 26321.51 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> L2 </TD> <TD align="right"> 16631.54 </TD> <TD align="right"> 0.91 </TD> <TD align="right"> 27769.46 </TD> <TD align="right"> 1.51 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 26918.98 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> fixed </TD> <TD align="right"> 11422.59 </TD> <TD align="right"> 0.62 </TD> <TD align="right"> 27028.94 </TD> <TD align="right"> 1.47 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 22992.45 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> L1 </TD> <TD align="right"> 15410.81 </TD> <TD align="right"> 0.84 </TD> <TD align="right"> 27970.24 </TD> <TD align="right"> 1.52 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 25154.28 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> L2 </TD> <TD align="right"> 15304.71 </TD> <TD align="right"> 0.83 </TD> <TD align="right"> 27395.68 </TD> <TD align="right"> 1.49 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 26559.39 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> fixed </TD> <TD align="right"> 4216.09 </TD> <TD align="right"> 0.23 </TD> <TD align="right"> 10796.44 </TD> <TD align="right"> 0.59 </TD> <TD align="right"> 18376.67 </TD> <TD align="right"> 9680.51 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> L1 </TD> <TD align="right"> 17158.98 </TD> <TD align="right"> 0.93 </TD> <TD align="right"> 28441.79 </TD> <TD align="right"> 1.54 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 27334.26 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> L2 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 28531.03 </TD> <TD align="right"> 1.54 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 28531.03 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> fixed </TD> <TD align="right"> 14886.65 </TD> <TD align="right"> 0.80 </TD> <TD align="right"> 27808.39 </TD> <TD align="right"> 1.50 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 25673.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> L1 </TD> <TD align="right"> 16502.83 </TD> <TD align="right"> 0.89 </TD> <TD align="right"> 28413.48 </TD> <TD align="right"> 1.54 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 26801.96 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> L2 </TD> <TD align="right"> 16701.71 </TD> <TD align="right"> 0.90 </TD> <TD align="right"> 28176.46 </TD> <TD align="right"> 1.52 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 27451.42 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> fixed </TD> <TD align="right"> 10456.70 </TD> <TD align="right"> 0.57 </TD> <TD align="right"> 27109.79 </TD> <TD align="right"> 1.47 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 23143.82 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> L1 </TD> <TD align="right"> 15191.71 </TD> <TD align="right"> 0.82 </TD> <TD align="right"> 28339.77 </TD> <TD align="right"> 1.53 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 25753.49 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> L2 </TD> <TD align="right"> 15193.26 </TD> <TD align="right"> 0.82 </TD> <TD align="right"> 27861.11 </TD> <TD align="right"> 1.51 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 26966.90 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> fixed </TD> <TD align="right"> 4334.74 </TD> <TD align="right"> 0.23 </TD> <TD align="right"> 10782.66 </TD> <TD align="right"> 0.58 </TD> <TD align="right"> 18495.32 </TD> <TD align="right"> 9666.74 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> L1 </TD> <TD align="right"> 17659.45 </TD> <TD align="right"> 0.92 </TD> <TD align="right"> 29754.66 </TD> <TD align="right"> 1.54 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 28390.30 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> L2 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 29884.61 </TD> <TD align="right"> 1.55 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 29884.61 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> fixed </TD> <TD align="right"> 15760.53 </TD> <TD align="right"> 0.82 </TD> <TD align="right"> 29271.16 </TD> <TD align="right"> 1.52 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 27189.40 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> L1 </TD> <TD align="right"> 16871.62 </TD> <TD align="right"> 0.88 </TD> <TD align="right"> 29675.05 </TD> <TD align="right"> 1.54 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 27770.16 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> L2 </TD> <TD align="right"> 17196.64 </TD> <TD align="right"> 0.89 </TD> <TD align="right"> 29469.12 </TD> <TD align="right"> 1.53 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 28395.76 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> fixed </TD> <TD align="right"> 11465.52 </TD> <TD align="right"> 0.59 </TD> <TD align="right"> 28761.45 </TD> <TD align="right"> 1.49 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 24711.55 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> L1 </TD> <TD align="right"> 15312.47 </TD> <TD align="right"> 0.79 </TD> <TD align="right"> 29516.09 </TD> <TD align="right"> 1.53 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 26520.45 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> L2 </TD> <TD align="right"> 15477.44 </TD> <TD align="right"> 0.80 </TD> <TD align="right"> 29060.32 </TD> <TD align="right"> 1.51 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 27708.94 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> fixed </TD> <TD align="right"> 6760.78 </TD> <TD align="right"> 0.35 </TD> <TD align="right"> 12609.11 </TD> <TD align="right"> 0.65 </TD> <TD align="right"> 19279.75 </TD> <TD align="right"> 10927.01 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
   </TABLE>

```r
writeLines(print(xtable(df)), "output.md")
```

% latex table generated in R 3.1.0 by xtable 1.7-3 package
% Tue Jun  3 15:05:26 2014
\begin{table}[ht]
\centering
\begin{tabular}{rlrrrrrrrrrr}
  \hline
 & penalty\_fn & ignore\_cost & ignore\_fraction & assume\_cost & assume\_fraction & normalize\_optimal\_free & normalize\_optimal\_cost & sigma\_g & reduction & L2 & L1 \\ 
  \hline
1 & L1 & 17187.64 & 0.94 & 28052.13 & 1.53 & 18376.67 & 26898.38 & 0.01 & 0.10 &   1 &   1 \\ 
  2 & L2 & 18376.67 & 1.00 & 28095.87 & 1.53 & 18376.67 & 28095.87 & 0.01 & 0.10 &   1 &   1 \\ 
  3 & fixed & 15284.62 & 0.83 & 27760.96 & 1.51 & 18376.67 & 25900.76 & 0.01 & 0.10 &   1 &   1 \\ 
  4 & L1 & 16593.12 & 0.90 & 28052.13 & 1.53 & 18376.67 & 26321.51 & 0.01 & 0.20 &   2 &   1 \\ 
  5 & L2 & 16631.54 & 0.91 & 27769.46 & 1.51 & 18376.67 & 26918.98 & 0.01 & 0.20 &   2 &   1 \\ 
  6 & fixed & 11422.59 & 0.62 & 27028.94 & 1.47 & 18376.67 & 22992.45 & 0.01 & 0.20 &   2 &   1 \\ 
  7 & L1 & 15410.81 & 0.84 & 27970.24 & 1.52 & 18376.67 & 25154.28 & 0.01 & 0.30 &   3 &   1 \\ 
  8 & L2 & 15304.71 & 0.83 & 27395.68 & 1.49 & 18376.67 & 26559.39 & 0.01 & 0.30 &   3 &   1 \\ 
  9 & fixed & 4216.09 & 0.23 & 10796.44 & 0.59 & 18376.67 & 9680.51 & 0.01 & 0.30 &   3 &   1 \\ 
  10 & L1 & 17158.98 & 0.93 & 28441.79 & 1.54 & 18495.32 & 27334.26 & 0.10 & 0.10 &   1 &   2 \\ 
  11 & L2 & 18495.32 & 1.00 & 28531.03 & 1.54 & 18495.32 & 28531.03 & 0.10 & 0.10 &   1 &   2 \\ 
  12 & fixed & 14886.65 & 0.80 & 27808.39 & 1.50 & 18495.32 & 25673.20 & 0.10 & 0.10 &   1 &   2 \\ 
  13 & L1 & 16502.83 & 0.89 & 28413.48 & 1.54 & 18495.32 & 26801.96 & 0.10 & 0.20 &   2 &   2 \\ 
  14 & L2 & 16701.71 & 0.90 & 28176.46 & 1.52 & 18495.32 & 27451.42 & 0.10 & 0.20 &   2 &   2 \\ 
  15 & fixed & 10456.70 & 0.57 & 27109.79 & 1.47 & 18495.32 & 23143.82 & 0.10 & 0.20 &   2 &   2 \\ 
  16 & L1 & 15191.71 & 0.82 & 28339.77 & 1.53 & 18495.32 & 25753.49 & 0.10 & 0.30 &   3 &   2 \\ 
  17 & L2 & 15193.26 & 0.82 & 27861.11 & 1.51 & 18495.32 & 26966.90 & 0.10 & 0.30 &   3 &   2 \\ 
  18 & fixed & 4334.74 & 0.23 & 10782.66 & 0.58 & 18495.32 & 9666.74 & 0.10 & 0.30 &   3 &   2 \\ 
  19 & L1 & 17659.45 & 0.92 & 29754.66 & 1.54 & 19279.75 & 28390.30 & 0.20 & 0.10 &   1 &   3 \\ 
  20 & L2 & 19279.75 & 1.00 & 29884.61 & 1.55 & 19279.75 & 29884.61 & 0.20 & 0.10 &   1 &   3 \\ 
  21 & fixed & 15760.53 & 0.82 & 29271.16 & 1.52 & 19279.75 & 27189.40 & 0.20 & 0.10 &   1 &   3 \\ 
  22 & L1 & 16871.62 & 0.88 & 29675.05 & 1.54 & 19279.75 & 27770.16 & 0.20 & 0.20 &   2 &   3 \\ 
  23 & L2 & 17196.64 & 0.89 & 29469.12 & 1.53 & 19279.75 & 28395.76 & 0.20 & 0.20 &   2 &   3 \\ 
  24 & fixed & 11465.52 & 0.59 & 28761.45 & 1.49 & 19279.75 & 24711.55 & 0.20 & 0.20 &   2 &   3 \\ 
  25 & L1 & 15312.47 & 0.79 & 29516.09 & 1.53 & 19279.75 & 26520.45 & 0.20 & 0.30 &   3 &   3 \\ 
  26 & L2 & 15477.44 & 0.80 & 29060.32 & 1.51 & 19279.75 & 27708.94 & 0.20 & 0.30 &   3 &   3 \\ 
  27 & fixed & 6760.78 & 0.35 & 12609.11 & 0.65 & 19279.75 & 10927.01 & 0.20 & 0.30 &   3 &   3 \\ 
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
## Time difference of 3.777 mins
```
