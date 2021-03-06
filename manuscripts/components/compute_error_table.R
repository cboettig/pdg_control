compute_error_table <- function(r){
  
  index <- dt_npv[,closest(r, (npv0-value)/npv0), by=variable]
  apples_index <- index$V1
  names(apples_index) <- index$variable
  
  L2_policy <- policies$L2[[apples_index["L2"]]]$D
  L1_policy <- policies$L1[[apples_index["L1"]]]$D
  fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D

  reps <- 1:500
  names(reps) = paste("rep", 1:length(reps), sep="_") # treat as a factor
  seeds <- 1:500
  sims <- list(
    L1 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                                                 L1_policy, z_g, z_m, z_i, 
                                                 opt$D, profit=profit, penalty=L1(c2[apples_index["L1"]]), seed=seeds[x])), 
    L2 = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                                                 L2_policy, z_g, z_m, z_i, 
                                                 opt$D, profit=profit, penalty=L2(c2[apples_index["L2"]]), seed=seeds[x])),
    fixed = lapply(reps, function(x) simulate_optim(f, pars, x_grid, h_grid, x0, 
                                                    fixed_policy, z_g, z_m, z_i, 
                                                    opt$D, profit=profit, penalty=fixed(c2[apples_index["fixed"]]), seed=seeds[x]))
    
  )


#Make data tidy (melt), fast (data.tables), and nicely labeled.
  dat <- melt(sims, id=names(sims[[1]][[1]]))  
  dt <- data.table(dat)
  setnames(dt, "L2", "replicate") # names are nice
  setnames(dt, "L1", "penalty_fn") # names are nice
  
  
  
  optimal_cost <- dt[, sum(profit_fishing - policy_cost), by=penalty_fn ] 
  ignore_when_present <- dt[, sum(profit_fishing_alt - policy_cost_alt), by=penalty_fn] 
  assume_when_absent <- dt[, sum(profit_fishing), by=penalty_fn]
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
  
  error_costs <- cbind(error_costs, sigma_g = sigma_g, reduction = r)
} 
 