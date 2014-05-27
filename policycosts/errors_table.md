

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
compute_error_table <- function(r = 2, sigma = 0.2){
reduction_range = c(.1, .2, .3)
reduction <- reduction_range[r]

price = 10
c0 = 0
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

sigma_g <- 0.3
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)

tmp <- lapply(1:100, function(i){ 
  set.seed(i)
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i)
  })
tmp <- melt(tmp, id = names(tmp[[1]]))
sum(tmp$profit)

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
## 1          L1       13284          0.7963       16678          0.9997
## 2          L2       16683          1.0000       16683          1.0000
## 3       fixed       12778          0.7660       16771          1.0053
## 4          L1       11988          0.7186       16472          0.9874
## 5          L2       11458          0.6868       15089          0.9045
## 6       fixed        8101          0.4856       16771          1.0053
## 7          L1        9062          0.5432       16148          0.9679
## 8          L2        6964          0.4174       14310          0.8577
## 9       fixed        3689          0.2211       15971          0.9573
## 10         L1       13744          0.7962       17226          0.9980
## 11         L2       17261          1.0000       17261          1.0000
## 12      fixed       14743          0.8541       17139          0.9929
## 13         L1       12080          0.6998       17213          0.9972
## 14         L2       11613          0.6728       15851          0.9183
## 15      fixed       11966          0.6933       16852          0.9763
## 16         L1        8864          0.5135       17207          0.9969
## 17         L2        7212          0.4178       14794          0.8571
## 18      fixed        9223          0.5343       16459          0.9535
## 19         L1       14973          0.8097       18451          0.9978
## 20         L2       18491          1.0000       18491          1.0000
## 21      fixed       14920          0.8068       18380          0.9940
## 22         L1       13242          0.7161       18360          0.9929
## 23         L2       12541          0.6782       17631          0.9535
## 24      fixed       11117          0.6012       18287          0.9889
## 25         L1        9777          0.5287       18237          0.9862
## 26         L2        8155          0.4410       16459          0.8901
## 27      fixed        6891          0.3727       18024          0.9747
##    normalize_optimal_free normalize_optimal_cost sigma_g reduction L2 L1
## 1                   16683                  13247    0.01       0.1  1  1
## 2                   16683                  16683    0.01       0.1  1  1
## 3                   16683                  14508    0.01       0.1  1  1
## 4                   16683                  11884    0.01       0.2  2  1
## 5                   16683                  13727    0.01       0.2  2  1
## 6                   16683                  11680    0.01       0.2  2  1
## 7                   16683                  10085    0.01       0.3  3  1
## 8                   16683                  13036    0.01       0.3  3  1
## 9                   16683                  10113    0.01       0.3  3  1
## 10                  17261                  13617    0.10       0.1  1  2
## 11                  17261                  17261    0.10       0.1  1  2
## 12                  17261                  14818    0.10       0.1  1  2
## 13                  17261                  11827    0.10       0.2  2  2
## 14                  17261                  14062    0.10       0.2  2  2
## 15                  17261                  12270    0.10       0.2  2  2
## 16                  17261                   8367    0.10       0.3  3  2
## 17                  17261                  13502    0.10       0.3  3  2
## 18                  17261                   9850    0.10       0.3  3  2
## 19                  18491                  15144    0.20       0.1  1  3
## 20                  18491                  18491    0.20       0.1  1  3
## 21                  18491                  15778    0.20       0.1  1  3
## 22                  18491                  13581    0.20       0.2  2  3
## 23                  18491                  14973    0.20       0.2  2  3
## 24                  18491                  13130    0.20       0.2  2  3
## 25                  18491                  10684    0.20       0.3  3  3
## 26                  18491                  14111    0.20       0.3  3  3
## 27                  18491                  10443    0.20       0.3  3  3
```



```r
#df <- df[1:5]
library(xtable)
print(xtable(df), type="html")
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Thu Nov  7 20:51:21 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> penalty_fn </TH> <TH> ignore_cost </TH> <TH> ignore_fraction </TH> <TH> assume_cost </TH> <TH> assume_fraction </TH> <TH> normalize_optimal_free </TH> <TH> normalize_optimal_cost </TH> <TH> sigma_g </TH> <TH> reduction </TH> <TH> L2 </TH> <TH> L1 </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> L1 </TD> <TD align="right"> 13381.38 </TD> <TD align="right"> 0.79 </TD> <TD align="right"> 16939.24 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 13354.13 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> L2 </TD> <TD align="right"> 11473.10 </TD> <TD align="right"> 0.68 </TD> <TD align="right"> 15368.64 </TD> <TD align="right"> 0.91 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 13833.72 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> fixed </TD> <TD align="right"> 14135.12 </TD> <TD align="right"> 0.84 </TD> <TD align="right"> 16875.95 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 14326.86 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> L1 </TD> <TD align="right"> 10312.53 </TD> <TD align="right"> 0.61 </TD> <TD align="right"> 16845.02 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 9843.99 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> L2 </TD> <TD align="right"> 7137.32 </TD> <TD align="right"> 0.42 </TD> <TD align="right"> 14433.86 </TD> <TD align="right"> 0.85 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 13305.76 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> fixed </TD> <TD align="right"> 10953.31 </TD> <TD align="right"> 0.65 </TD> <TD align="right"> 16584.22 </TD> <TD align="right"> 0.98 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 11525.64 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> L1 </TD> <TD align="right"> 6179.21 </TD> <TD align="right"> 0.37 </TD> <TD align="right"> 16788.42 </TD> <TD align="right"> 0.99 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 4713.25 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> L2 </TD> <TD align="right"> -7441.61 </TD> <TD align="right"> -0.44 </TD> <TD align="right"> 13301.51 </TD> <TD align="right"> 0.79 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 12625.07 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> fixed </TD> <TD align="right"> 8034.52 </TD> <TD align="right"> 0.47 </TD> <TD align="right"> 16024.94 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 16916.94 </TD> <TD align="right"> 9327.57 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> L1 </TD> <TD align="right"> 14973.02 </TD> <TD align="right"> 0.81 </TD> <TD align="right"> 18451.17 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 15143.69 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> L2 </TD> <TD align="right"> 12541.29 </TD> <TD align="right"> 0.68 </TD> <TD align="right"> 17630.65 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 14973.26 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> fixed </TD> <TD align="right"> 14465.87 </TD> <TD align="right"> 0.78 </TD> <TD align="right"> 18420.03 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 15547.30 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> L1 </TD> <TD align="right"> 11526.93 </TD> <TD align="right"> 0.62 </TD> <TD align="right"> 18274.97 </TD> <TD align="right"> 0.99 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 12166.27 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> L2 </TD> <TD align="right"> 8154.80 </TD> <TD align="right"> 0.44 </TD> <TD align="right"> 16458.95 </TD> <TD align="right"> 0.89 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 14111.15 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> fixed </TD> <TD align="right"> 10297.38 </TD> <TD align="right"> 0.56 </TD> <TD align="right"> 18259.56 </TD> <TD align="right"> 0.99 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 12603.00 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> L1 </TD> <TD align="right"> 6440.60 </TD> <TD align="right"> 0.35 </TD> <TD align="right"> 17885.39 </TD> <TD align="right"> 0.97 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 8108.31 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> L2 </TD> <TD align="right"> -5229.90 </TD> <TD align="right"> -0.28 </TD> <TD align="right"> 14894.07 </TD> <TD align="right"> 0.81 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 13179.37 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> fixed </TD> <TD align="right"> 5781.83 </TD> <TD align="right"> 0.31 </TD> <TD align="right"> 17929.89 </TD> <TD align="right"> 0.97 </TD> <TD align="right"> 18491.32 </TD> <TD align="right"> 9616.35 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   2 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> L1 </TD> <TD align="right"> 19803.90 </TD> <TD align="right"> 0.83 </TD> <TD align="right"> 23743.28 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> 19906.24 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> L2 </TD> <TD align="right"> 14900.10 </TD> <TD align="right"> 0.63 </TD> <TD align="right"> 23551.58 </TD> <TD align="right"> 0.99 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> 17817.72 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> fixed </TD> <TD align="right"> 19636.82 </TD> <TD align="right"> 0.83 </TD> <TD align="right"> 23782.02 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> 20542.02 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.10 </TD> <TD align="right">   1 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> L1 </TD> <TD align="right"> 15938.52 </TD> <TD align="right"> 0.67 </TD> <TD align="right"> -299999976471.53 </TD> <TD align="right"> -12645374.51 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> -299999983816.73 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> L2 </TD> <TD align="right"> 8584.62 </TD> <TD align="right"> 0.36 </TD> <TD align="right"> 22496.08 </TD> <TD align="right"> 0.95 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> 16627.19 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> fixed </TD> <TD align="right"> 14665.50 </TD> <TD align="right"> 0.62 </TD> <TD align="right"> -299999976708.30 </TD> <TD align="right"> -12645374.51 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> -299999983164.87 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.20 </TD> <TD align="right">   2 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> L1 </TD> <TD align="right"> 10216.03 </TD> <TD align="right"> 0.43 </TD> <TD align="right"> -299999976767.67 </TD> <TD align="right"> -12645374.52 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> -299999988753.36 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> L2 </TD> <TD align="right"> -10584.43 </TD> <TD align="right"> -0.45 </TD> <TD align="right"> 19839.12 </TD> <TD align="right"> 0.84 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> 14962.22 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> fixed </TD> <TD align="right"> 9476.01 </TD> <TD align="right"> 0.40 </TD> <TD align="right"> -599999977068.36 </TD> <TD align="right"> -25290750.03 </TD> <TD align="right"> 23724.09 </TD> <TD align="right"> -599999986390.38 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.30 </TD> <TD align="right">   3 </TD> <TD align="right">   3 </TD> </TR>
   </TABLE>

```r
writeLines(print(xtable(df)), "output.md")
```

% latex table generated in R 3.0.2 by xtable 1.7-1 package
% Thu Nov  7 20:51:21 2013
\begin{table}[ht]
\centering
\begin{tabular}{rlrrrrrrrrrr}
  \hline
 & penalty\_fn & ignore\_cost & ignore\_fraction & assume\_cost & assume\_fraction & normalize\_optimal\_free & normalize\_optimal\_cost & sigma\_g & reduction & L2 & L1 \\ 
  \hline
1 & L1 & 13381.38 & 0.79 & 16939.24 & 1.00 & 16916.94 & 13354.13 & 0.05 & 0.10 &   1 &   1 \\ 
  2 & L2 & 11473.10 & 0.68 & 15368.64 & 0.91 & 16916.94 & 13833.72 & 0.05 & 0.10 &   1 &   1 \\ 
  3 & fixed & 14135.12 & 0.84 & 16875.95 & 1.00 & 16916.94 & 14326.86 & 0.05 & 0.10 &   1 &   1 \\ 
  4 & L1 & 10312.53 & 0.61 & 16845.02 & 1.00 & 16916.94 & 9843.99 & 0.05 & 0.20 &   2 &   1 \\ 
  5 & L2 & 7137.32 & 0.42 & 14433.86 & 0.85 & 16916.94 & 13305.76 & 0.05 & 0.20 &   2 &   1 \\ 
  6 & fixed & 10953.31 & 0.65 & 16584.22 & 0.98 & 16916.94 & 11525.64 & 0.05 & 0.20 &   2 &   1 \\ 
  7 & L1 & 6179.21 & 0.37 & 16788.42 & 0.99 & 16916.94 & 4713.25 & 0.05 & 0.30 &   3 &   1 \\ 
  8 & L2 & -7441.61 & -0.44 & 13301.51 & 0.79 & 16916.94 & 12625.07 & 0.05 & 0.30 &   3 &   1 \\ 
  9 & fixed & 8034.52 & 0.47 & 16024.94 & 0.95 & 16916.94 & 9327.57 & 0.05 & 0.30 &   3 &   1 \\ 
  10 & L1 & 14973.02 & 0.81 & 18451.17 & 1.00 & 18491.32 & 15143.69 & 0.20 & 0.10 &   1 &   2 \\ 
  11 & L2 & 12541.29 & 0.68 & 17630.65 & 0.95 & 18491.32 & 14973.26 & 0.20 & 0.10 &   1 &   2 \\ 
  12 & fixed & 14465.87 & 0.78 & 18420.03 & 1.00 & 18491.32 & 15547.30 & 0.20 & 0.10 &   1 &   2 \\ 
  13 & L1 & 11526.93 & 0.62 & 18274.97 & 0.99 & 18491.32 & 12166.27 & 0.20 & 0.20 &   2 &   2 \\ 
  14 & L2 & 8154.80 & 0.44 & 16458.95 & 0.89 & 18491.32 & 14111.15 & 0.20 & 0.20 &   2 &   2 \\ 
  15 & fixed & 10297.38 & 0.56 & 18259.56 & 0.99 & 18491.32 & 12603.00 & 0.20 & 0.20 &   2 &   2 \\ 
  16 & L1 & 6440.60 & 0.35 & 17885.39 & 0.97 & 18491.32 & 8108.31 & 0.20 & 0.30 &   3 &   2 \\ 
  17 & L2 & -5229.90 & -0.28 & 14894.07 & 0.81 & 18491.32 & 13179.37 & 0.20 & 0.30 &   3 &   2 \\ 
  18 & fixed & 5781.83 & 0.31 & 17929.89 & 0.97 & 18491.32 & 9616.35 & 0.20 & 0.30 &   3 &   2 \\ 
  19 & L1 & 19803.90 & 0.83 & 23743.28 & 1.00 & 23724.09 & 19906.24 & 0.50 & 0.10 &   1 &   3 \\ 
  20 & L2 & 14900.10 & 0.63 & 23551.58 & 0.99 & 23724.09 & 17817.72 & 0.50 & 0.10 &   1 &   3 \\ 
  21 & fixed & 19636.82 & 0.83 & 23782.02 & 1.00 & 23724.09 & 20542.02 & 0.50 & 0.10 &   1 &   3 \\ 
  22 & L1 & 15938.52 & 0.67 & -299999976471.53 & -12645374.51 & 23724.09 & -299999983816.73 & 0.50 & 0.20 &   2 &   3 \\ 
  23 & L2 & 8584.62 & 0.36 & 22496.08 & 0.95 & 23724.09 & 16627.19 & 0.50 & 0.20 &   2 &   3 \\ 
  24 & fixed & 14665.50 & 0.62 & -299999976708.30 & -12645374.51 & 23724.09 & -299999983164.87 & 0.50 & 0.20 &   2 &   3 \\ 
  25 & L1 & 10216.03 & 0.43 & -299999976767.67 & -12645374.52 & 23724.09 & -299999988753.36 & 0.50 & 0.30 &   3 &   3 \\ 
  26 & L2 & -10584.43 & -0.45 & 19839.12 & 0.84 & 23724.09 & 14962.22 & 0.50 & 0.30 &   3 &   3 \\ 
  27 & fixed & 9476.01 & 0.40 & -599999977068.36 & -25290750.03 & 23724.09 & -599999986390.38 & 0.50 & 0.30 &   3 &   3 \\ 
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
## Time difference of 2.22 mins
```

