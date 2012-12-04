
fig4 <- function(fraction_lost, ...){
closest <- function(x, v){
  which.min(abs(v-x))
}
dt_npv <- data.table(fees)
index <- dt_npv[,closest(fraction_lost, (npv0-value)/npv0), by=variable]
apples_index <- index$V1
names(apples_index) = index$variable
apples <- c2[index$V1]
names(apples) = index$variable
apples

L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D
free_increase_policy <- policies$free_increase[[apples_index["free_increase"]]]$D
free_decrease_policy <- policies$free_decrease[[apples_index["free_decrease"]]]$D
quad_policy <- policies$quad[[apples_index["quad"]]]$D

quad_profit <- profit_harvest(price = price, c0 = c0, c1 = apples["quad"]) 
sims <- list(
  L1 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L1_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L1(apples["L1"]), seed=seed), 
  L2 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L2_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L2(apples["L2"]), seed=seed),
  fixed = simulate_optim(f, pars, x_grid, h_grid, x0, 
                         fixed_policy, z_g, z_m, z_i, 
                         opt$D, profit=profit, penalty=fixed(apples["fixed"]), seed=seed),
  increase = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            free_increase_policy, z_g, z_m, z_i, 
                            opt$D, profit=profit, penalty= free_increase(apples["increase"]), seed=seed),
  decrease = simulate_optim(f, pars, x_grid, h_girpdgd, x0, 
                            free_decrease_policy, z_g, z_m, z_i, 
                            opt$D, profit=profit, penalty= free_decrease(apples["decrease"]), seed=seed),
  quad = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            quad_policy, z_g, z_m, z_i, 
                            opt$D, profit=quad_profit, penalty= none, seed=seed)
)

#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "penalty_fn") # names are nice


v <- dt[,var(harvest), by="penalty_fn"]
var <- v$V1
names(var) <- v$penalty_fn
acor <- dt[,acf(harvest, plot=F)$acf[2], by="penalty_fn"]$V1
names(acor) <- names(var)
out <- rbind(variance=var, autocorrelation=acor)
out
}