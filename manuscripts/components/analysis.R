
#### SETUP #######

parallel = TRUE
ncpu = 16


## ----profit_model--------------------------------------------------------
price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)

## ----c2_grid-------------------------------------------------------------
c2 <- exp(seq(0, log(41), length.out = 40))-1
c2 <- seq(0, 40, length.out=100)


## ----reduction-----------------------------------------------------------
reduction <- 0.25 


## ----setup---------------------------------------------------------------
seed <- 123                 # Random seed (replicable results)
delta <- 0.05               # economic discounting rate
OptTime <- 50               # stopping time
gridsize <- 50              # grid size for fish stock and harvest rate (discretized population)
sigma_g <- 0.2              # Noise in population growth
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


## ----reed--------------------------
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)


## ----fees----------------------------------------------------------------
L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
fixed <-  function(c2) function(h, h_prev) c2 * as.numeric( !(h == h_prev) )
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2
none <- function(h, h_prev)  0
penaltyfns <- list(L2=L2, L1=L1, fixed=fixed)




## ----parallel, include=FALSE---------------------------------------------
sfInit(cpu=ncpu, parallel=parallel)
sfLibrary(pdgControl)
sfExportAll()

##### ANALYSIS #############


## Apples to Apples comparison (Figure 2)
policies <- lapply(penaltyfns, function(penalty){
  sfLapply(c2, function(c2){
    policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                           profit, delta, reward, penalty = penalty(c2))
  }
  )
})

i <- which(x_grid > K)[1]
fees <- 
  lapply(policies, function(penalty) 
    sapply(penalty, function(c2_run)
      max(c2_run$V[i,]) # Would be penalty_free_V originally   ## this isn't correct for asym cases
    )
  )


## ----npv-plot-------------------------------------
npv0 <- max(fees$L1) # all have same max, at c2=0 
fees <- data.frame(c2=c2,fees)
fees <- melt(fees, id="c2")
fees <- subset(fees, variable %in% c("L1", "L2", "fixed"))

## ----apples-----------------------
closest <- function(x, v){
  which.min(abs(v-x))
}
dt_npv <- data.table(fees)
index <- dt_npv[,closest(reduction, (npv0-value)/npv0), by=variable]
apples_index <- index$V1
names(apples_index) <- index$variable
apples <- c2[index$V1]
names(apples) <- index$variable


## ----print_npv---------------------------------------
setkey(dt_npv, variable)
values <- apply(index, 1, function(x) dt_npv[x[1], ][as.integer(x[2]), ]$value)
percent.error <- (values - ((1-reduction)*npv0)) / ((1-reduction)*npv0)* 100
print_npv <- data.frame(model=index$variable, "value realized"=values, 
                        "percent of npv0" = 100*values/npv0, "percent error"=percent.error)



## ----Apply the apples-to-apples calibration------------------------------------
L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D


## ----Simulations for Figure 3-------------------
reps <- 1:500
names(reps) = paste("rep", 1:length(reps), sep="_") # treat as a factor
seeds <- 1:500
sims <- list(
  L1 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                                               L1_policy, z_g, z_m, z_i, 
                                               opt$D, profit=profit, penalty=L1(apples["L1"]), seed=seeds[x])), 
  L2 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                                               L2_policy, z_g, z_m, z_i, 
                                               opt$D, profit=profit, penalty=L2(apples["L2"]), seed=seeds[x])),
  fixed = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                                                  fixed_policy, z_g, z_m, z_i, 
                                                  opt$D, profit=profit, penalty=fixed(apples["fixed"]), seed=seeds[x]))
  
)


## ----tidy, dependson="simulate_policy"-----------------------------------
#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]][[1]]))  
dt <- data.table(dat)
setnames(dt, "L2", "replicate") # names are nice
setnames(dt, "L1", "penalty_fn") # names are nice




#####------------------Simulations for Figure 5---------------------------
## ----profit_calcs, dependson="tidy"--------------------------------------
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
error_costs <- merge(ignore_when_present, assume_when_absent, "penalty_fn")

# print_npv                             # theoretically acheivable profits under these costs
# optimal_cost$V1/optimal_free$V1       # actually realized 

error_costs <- cbind(error_costs, sigma_g = sigma_g, reduction = reduction)



## ----helper_fn_1---------------------------------------------------------
fig4 <- function(fraction_lost){
  closest <- function(x, v){
    which.min(abs(v-x))
  }
  dt_npv <- data.table(fees)
  index <- dt_npv[,closest(fraction_lost, (npv0-value)/npv0), by=variable]
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
  sims <- lapply(1:500, function(reps) list(
    L1 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                        L1_policy, z_g, z_m, z_i, 
                        opt$D, profit=profit, penalty=L1(apples["L1"])), 
    L2 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                        L2_policy, z_g, z_m, z_i, 
                        opt$D, profit=profit, penalty=L2(apples["L2"])),
    fixed = simulate_optim(f, pars, x_grid, h_grid, x0, 
                           fixed_policy, z_g, z_m, z_i, 
                           opt$D, profit=profit, penalty=fixed(apples["fixed"]))
    #  increase = simulate_optim(f, pars, x_grid, h_grid, x0, 
    #                            free_increase_policy, z_g, z_m, z_i, 
    #                            opt$D, profit=profit, penalty= free_increase(apples["increase"])),
    #  decrease = simulate_optim(f, pars, x_grid, h_grid, x0, 
    #                            free_decrease_policy, z_g, z_m, z_i, 
    #                            opt$D, profit=profit, penalty= free_decrease(apples["decrease"])),
    #  quad = simulate_optim(f, pars, x_grid, h_grid, x0, 
    #                            quad_policy, z_g, z_m, z_i, 
    #                            opt$D, profit=quad_profit, penalty= none)
  ))
  
  #Make data tidy (melt), fast (data.tables), and nicely labeled.
  dat <- melt(sims, id=names(sims[[1]][[1]]))  
  dt <- data.table(dat)
  setnames(dt, "L1", "reps") # names are nice
  setnames(dt, "L2", "penalty_fn") # names are nice
  
  dt
}


##### Figure 4 #########

## ----sim_at_each_apple---------------------------------------------------
frac_lost <- seq(0,1, length=20)
sims_at_each_apple <- lapply(frac_lost, fig4)


## ----stats-defs----------------------------------------------------------
summary_stats <- function(dt){
  harvest_var <- dt[,var(harvest), by=c("penalty_fn", "reps")]$V1
  harvest_acor <- dt[,acf(harvest, plot=F)$acf[2], by=c("penalty_fn", "reps")]$V1
  stock_var <- dt[,var(fishstock), by=c("penalty_fn", "reps")]$V1
  stock_acor <- dt[,acf(fishstock, plot=F)$acf[2], by=c("penalty_fn", "reps")]$V1
  xcor <- dt[, ccf(fishstock, harvest, plot=FALSE)$acf[2], by=c("penalty_fn", "reps")]
  setnames(xcor, "V1", "cross-correlation")
  
  df <- data.frame('stock.variance' = stock_var, 
                   'stock.autocorrelation' = stock_acor, 
                   'harvest.variance' = harvest_var, 
                   'harvest.autocorrelation' = harvest_acor, 
                   xcor)
  
  # long format
  melt(df, id=c("penalty_fn", "reps"))
}


## ----summary-stat-calc, dependson="sim_at_each_apple"--------------------
stats_summaries <- lapply(sims_at_each_apple, summary_stats)
# reformat list of data-frames as data frame with apple coef as a factor
stats_df <- melt(stats_summaries, id=names(stats_summaries[[1]]))
stats_df <- as_data_frame(stats_df)
stats_df$L1 <- frac_lost[stats_df$L1]
stats_df %>% 
  dplyr::rename(penalty_fraction = L1) -> stats_df


stats_df %>% 
  filter(variable != "cross.correlation") %>%
  separate(variable, c("measurement", "statistic"), sep="\\.") %>%
  filter(measurement == 'harvest') %>% 
  spread(statistic, value) ->
  f4df


## ----histograms----------------------------------------------------------
profits <- dt[, sum(profit_fishing), by=c('penalty_fn', 'replicate') ]
costs <- dt[, sum(policy_cost), by=c('penalty_fn', 'replicate') ]
reed_profits <- dt[, sum(profit_fishing_alt), by=c('penalty_fn', 'replicate') ]
reed_costs <- dt[, sum(policy_cost_alt), by=c('penalty_fn', 'replicate') ]
setnames(profits, "V1", "profits")
hist_dat <- melt(cbind(profits, costs = costs$V1, 
                       reed_profits = reed_profits$V1, reed_costs = reed_costs$V1),
                 id = c("penalty_fn", "replicate"))

## ----errors-table--------------------------------------------------------

source("components/compute_error_table.R")

reduction_list <- c(.1, .15, .2, .25, .3)
out <- lapply(reduction_list, function(r)
  compute_error_table(r = r))   
who <- names(out[[1]])
error_df <- melt(out, id=who)       


## ----mismatches--------------------------------------------------------

source("components/compute_mismatches.R")

reduction_list <- c(.1, .15, .2, .25)
out <- lapply(reduction_list, function(r)
  compute_mismatches(r = r))   
who <- names(out[[1]])
mismatches_df <- melt(out, id=who)       

