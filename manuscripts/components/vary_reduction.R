

robust_reduction <- function(reduction){

i <- which(x_grid > K)[1]
fees <- 
  lapply(policies, function(penalty) 
    sapply(penalty, function(c2_run)
      max(c2_run$V[i,]) # Would be penalty_free_V originally   ## this isn't correct for asym cases
    )
  )


## ----npv--------------------------------------
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

## ----Apply the apples-to-apples calibration------------------------------------
L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D


## ----Simulations for Figure 3-------------------
reps <- 1:500
names(reps) = paste("rep", 1:length(reps), sep="_") # treat as a factor
seeds <- 1:length(reps)
sims <- list(
  L1 = lapply(reps, 
              function(x) 
              simulate_optim(f, pars, x_grid, h_grid, x0, 
                             L1_policy, z_g, z_m, z_i, 
                             opt$D, profit=profit, 
                             penalty=L1(apples["L1"]), 
                             seed=seeds[x])), 
  L2 = lapply(reps, 
              function(x) 
              simulate_optim(f, pars, x_grid, h_grid, x0, 
                             L2_policy, z_g, z_m, z_i, 
                             opt$D, profit=profit, 
                             penalty=L2(apples["L2"]), 
                             seed=seeds[x])),
  fixed = lapply(reps, 
                 function(x) 
                 simulate_optim(f, pars, x_grid, h_grid, x0, 
                                fixed_policy, z_g, z_m, z_i, 
                                opt$D, profit=profit, 
                                penalty=fixed(apples["fixed"]), 
                                seed=seeds[x]))
  
)


## ----tidy, dependson="simulate_policy"-----------------------------------
#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]][[1]]))  
dt <- data.table(dat)
setnames(dt, "L2", "replicate") # names are nice
setnames(dt, "L1", "penalty_fn") # names are nice





##### Figure 4 #########
## ----sim_at_each_apple---------------------------------------------------
frac_lost <- seq(0,1, length=40)
sims_at_each_apple <- lapply(frac_lost, fig4)

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



### Figure 3
labeller <- function(variable,value){
    return(relabel[paste(value)])
}

dt3 <- dt %>% 
  filter(replicate=="rep_17" & time < OptTime) %>%
  select(time, harvest, harvest_alt, penalty_fn) %>%
  gather(variable, value, -time, -penalty_fn) 

dt3$penalty_fn <- factor(dt3$penalty_fn, levels=c("L1", "L2", "fixed"))

profits <- dt[, sum(profit_fishing), by=c('penalty_fn', 'replicate') ]
costs <- dt[, sum(policy_cost), by=c('penalty_fn', 'replicate') ]
reed_profits <- dt[, sum(profit_fishing_alt), by=c('penalty_fn', 'replicate') ]
reed_costs <- dt[, sum(policy_cost_alt), by=c('penalty_fn', 'replicate') ]
setnames(profits, "V1", "profits")
setnames(reed_profits, "V1", "profits")

Reed <- cbind(reed_profits, costs = reed_costs$V1, Assumption = "No adjustment penalty") 
Adj <- cbind(profits, costs = costs$V1, Assumption = "Adjustment penalty")

hist_dat <- melt(rbind(Adj, Reed), id=c("penalty_fn", "replicate", "Assumption"))

list(dt3 = dt3, dt5 = hist_dat) 

}


low <- robust_reduction(0.10)
high <- robust_reduction(0.30)


