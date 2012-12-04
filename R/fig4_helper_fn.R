
#' @export
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

L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D
free_increase_policy <- policies$free_increase[[apples_index["free_increase"]]]$D
free_decrease_policy <- policies$free_decrease[[apples_index["free_decrease"]]]$D
quad_policy <- policies$quad[[apples_index["quad"]]]$D

quad_profit <- profit_harvest(price = price, c0 = c0, c1 = apples["quad"]) 
sims <- lapply(1:50, function(reps) list(
  L1 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L1_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L1(apples["L1"])), 
  L2 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L2_policy, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L2(apples["L2"])),
  fixed = simulate_optim(f, pars, x_grid, h_grid, x0, 
                         fixed_policy, z_g, z_m, z_i, 
                         opt$D, profit=profit, penalty=fixed(apples["fixed"])),
  increase = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            free_increase_policy, z_g, z_m, z_i, 
                            opt$D, profit=profit, penalty= free_increase(apples["increase"])),
  decrease = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            free_decrease_policy, z_g, z_m, z_i, 
                            opt$D, profit=profit, penalty= free_decrease(apples["decrease"])),
  quad = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            quad_policy, z_g, z_m, z_i, 
                            opt$D, profit=quad_profit, penalty= none)
))

#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]][[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
setnames(dt, "L2", "penalty_fn") # names are nice

v <- dt[,var(harvest), by=c("penalty_fn", "reps")]
acor <- dt[,acf(harvest, plot=F)$acf[2], by=c("penalty_fn", "reps")]

df <- cbind(v, acor$V1)
setnames(df, c("V1", "V2"), c("var", "acor")) # names are nice


#Expected_var <- v[,mean(V1), by="penalty_fn"]
#Error_var <- v[,sd(V1), by="penalty_fn"]
#Expected_acor <- acor[,mean(V1), by="penalty_fn"]
#Error_acor <- acor[,sd(V1), by="penalty_fn"]
#out <- rbind(variance=var, autocorrelation=acor)
df
}