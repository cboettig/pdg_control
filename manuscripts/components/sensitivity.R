parallel = FALSE
ncpu = 4

## ----dependencies, include=FALSE-----------------------------------------
source("components/install.R")


## ----caching, include=FALSE----------------------------------------------
library("methods")
library("knitr")
basename <- "manuscript"
opts_chunk$set(fig.path = paste("components/figure/", basename, "-", sep=""),
               cache.path = paste("components/cache/", basename, "/", sep=""))
opts_chunk$set(cache = 2)
opts_chunk$set(tidy=FALSE, warning=FALSE, message=FALSE, 
               comment = NA, verbose = TRUE, echo=FALSE)

# PDF-based figures
opts_chunk$set(dev='pdf')


## ----libraries, include=FALSE, cache=FALSE-------------------------------
rm(list=ls())
library("rmarkdown")
library("pdgControl")
library("reshape2")
library("ggplot2")
library("data.table")
library("pander")
library("cboettigR")
library("ggthemes")
library("snowfall")
theme_set(theme_tufte())


## ----profit_model--------------------------------------------------------
price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 1)


## ----c2_grid-------------------------------------------------------------
c2 <- exp(seq(0, log(41), length.out = 40))-1
c2 <- seq(0, 40, length.out=100)


## ----reduction-----------------------------------------------------------
reduction <- 0.25 


## ----setup---------------------------------------------------------------
seed <- 123                 # Random seed (replicable results)
delta <- 0.05               # economic discounting rate
OptTime <- 20               # stopping time
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


## ----reed, dependson=c("setup", "profit_model")--------------------------
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


## ----bigloop, dependson=c("setup", "reed", "fees", "profit_model", "c2_grid")----
policies <- lapply(penaltyfns, function(penalty){
  sfLapply(c2, function(c2){
      policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                   profit, delta, reward, penalty = penalty(c2))
      }
  )
})


## ----, dependson="quadcosts"---------------------------------------------
i <- which(x_grid > K)[1]
fees <- 
lapply(policies, function(penalty) 
  sapply(penalty, function(c2_run)
    max(c2_run$V[i,]) # Would be penalty_free_V originally   ## this isn't correct for asym cases
  )
)


## ----npv-plot, dependson="quadcosts"-------------------------------------
npv0 <- max(fees$L1) # all have same max, at c2=0 
fees <- data.frame(c2=c2,fees)
fees <- melt(fees, id="c2")
fees <- subset(fees, variable %in% c("L1", "L2", "fixed"))

#ggplot(fees, aes(c2, value, col=variable)) + geom_point() + geom_line()


## ----apples_plot, dependson="quadcosts"----------------------------------
#ggplot(fees, aes(c2, (npv0-value)/npv0, col=variable)) + geom_point() + geom_line()


## ----apples, dependson=c("quadcosts", "reduction")-----------------------
closest <- function(x, v){
  which.min(abs(v-x))
}
dt_npv <- data.table(fees)
index <- dt_npv[,closest(reduction, (npv0-value)/npv0), by=variable]
apples_index <- index$V1
names(apples_index) = index$variable
apples <- c2[index$V1]


## ----reductions----------------------------------------------------------
reduction_table <- sapply(c("0.10" = 0.10, "0.20" = 0.20, "0.30" = 0.30), function(reduction){
index <- dt_npv[,closest(reduction, (npv0-value)/npv0), by=variable]
apples <- c2[index$V1]
names(apples) = index$variable
apples })
write.csv(reduction_table, "../data/reduction_table.csv")


## ----print_npv, dependson="apples"---------------------------------------
setkey(dt_npv, variable)
values <- apply(index, 1, function(x) dt_npv[x[1], ][as.integer(x[2]), ]$value)
percent.error <- (values - ((1-reduction)*npv0)) / ((1-reduction)*npv0)* 100
print_npv <- data.frame(model=index$variable, "value realized"=values, 
                        "percent of npv0" = 100*values/npv0, "percent error"=percent.error)




## ----npv_table, dependson="print_npv", results='asis'--------------------
# verify that we have fine enough c_2 sampling to hit about 75% NPV0 on each penalty form
# pandoc.table(print_npv)


## ----policynames---------------------------------------------------------
L2_policy <- policies$L2[[apples_index["L2"]]]$D
L1_policy <- policies$L1[[apples_index["L1"]]]$D
fixed_policy <- policies$fixed[[apples_index["fixed"]]]$D


## ----simulate_policy, dependson=c("policynames", "apples")---------------
reps <- 1:100
names(reps) = paste("rep", 1:length(reps), sep="_") # treat as a factor
seeds <- 1:100
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



## ------------------------------------------------------------------------
v <- dt[,var(harvest), by="penalty_fn"]
var <- v$V1
names(var) <- v$penalty_fn
acor <- dt[,acf(harvest, plot=F)$acf[2], by="penalty_fn"]$V1
names(acor) <- names(var)
out <- rbind(var=var, a=acor)


## ------------------------------------------------------------------------
frac_lost <- seq(0,1, length=20)


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
sims <- lapply(1:50, function(reps) list(
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


## ----helper_fn_2---------------------------------------------------------
# Helper functions to extract the summary stats in different variables
stats_harvest <- function(dt){
v <- dt[,var(harvest), by=c("penalty_fn", "reps")]
acor <- dt[,acf(harvest, plot=F)$acf[2], by=c("penalty_fn", "reps")]
df <- cbind(v, acor$V1)
setnames(df, c("V1", "V2"), c("var", "acor")) # names are nice
df
}

stats_fishstock <- function(dt){
v <- dt[,var(fishstock), by=c("penalty_fn", "reps")]
acor <- dt[,acf(fishstock, plot=F)$acf[2], by=c("penalty_fn", "reps")]
df <- cbind(v, acor$V1)
setnames(df, c("V1", "V2"), c("var", "acor")) # names are nice
df
}


## ----helper_fn_3---------------------------------------------------------
#' a simple function for reorganizing the data over the different "faction-lost" levels
get_trends <- function(tmp2){
  out <- melt(tmp2, id=c("reps", "var", "acor"))
  colnames(out) = c("reps", "var", "acor", "nothing", "penalty", "index")
  out <- cbind(out[c(1,2,3,5)], fraction = frac_lost[out$index])
  out <- data.table(out)
  Ev = out[,mean(var),by=c('penalty', 'fraction')]
  SDv = out[,sd(var),by=c('penalty', 'fraction')]
  Ea = out[,mean(acor),by=c('penalty', 'fraction')]
  SDa = out[,sd(acor),by=c('penalty', 'fraction')]
  harvest_trends <- data.table(penalty = Ev$penalty, 
                               fraction = Ev$fraction, 
                               Ev = Ev$V1, Ea = Ea$V1, 
                               SDv=SDv$V1, SDa = SDa$V1)
}


## ----Figure4_harvest-----------------------------------------------------
sims_at_each_apple <- lapply(frac_lost, fig4)
harvest_stats <- lapply(sims_at_each_apple, stats_harvest)
harvest_trends <- get_trends(harvest_stats)


## ----fishstock_trends----------------------------------------------------
fishstock_stats <- lapply(sims_at_each_apple, stats_fishstock)
fishstock_trends <- get_trends(fishstock_stats)


## ----Figure_1, fig.cap="Expected net present value by functional form of penalty. Horizontal axis shows the coefficient $c_i$ governing the magnitude of the policy cost, while vertical axis shows fraction of the maximum expected net present value dissipated by the policy cost, $\\tfrac{NPV_0 - NPV_i}{NPV_0}$.  The horizontal line indicates a value of the stock that is reduced by `r reduction * 100`% from the maximum expected value in the absence of policy adjustment costs, $NPV_0({\\bf h_0^*})$. Selecting the coefficient $c_i$ corresponding to this value in each functional form allows us to make consistent comparisons across the different functional forms of policy costs.", dependson="quadcosts", fig.width=4, fig.height=3----

ggplot(fees, aes(c2, (npv0-value)/npv0, lty=variable)) + 
  geom_line() + 
  geom_hline(aes(hintercept=reduction), linetype=4) + 
  xlab("penalty coefficient") +
  ylab("Reduction in net present value") + 
  theme_tufte() 



## ----Figure_2, dependson="tidy", fig.cap="Example realization of optimal harvesting strategy and resulting fish stock sizes under the different functional forms of adjustment costs.", fig.height=3, fig.width=8----

fig2_df <- as.data.frame(subset(dt, replicate=='rep_2'))
fig2_df <- fig2_df[c("time", "fishstock", "alternate", "harvest", "harvest_alt", "penalty_fn")]
fig2_df <- melt(fig2_df, id = c("time", "penalty_fn"))
fig2_df <- data.frame(fig2_df, baseline = fig2_df$variable)
variable_map <- c(fishstock = "fish_stock", alternate = "fish_stock", harvest = "harvest", harvest_alt = "harvest")
baseline_map <- c(fishstock = "penalty", alternate = "no_penalty", harvest = "penalty", harvest_alt = "no_penalty")
fig2_df$variable <- variable_map[fig2_df$variable]
fig2_df$baseline <- baseline_map[fig2_df$baseline]


ggplot(fig2_df, aes(time, value, col=baseline)) +
  geom_line() +
  facet_grid(variable~penalty_fn) + 
  labs(x="time", y="stock size", title = "Example Stock & Harvest Dynamics")  + 
  theme_tufte()




## ----Figure_4, fig.cap="(a) Variance of optimal harvest rates through time as a function of size of adjustment penalties, the latter being calculated as percentage of the maximum expected NPV in the basic case with no policy adjustment costs. (b) Autocorrelation of optimal harvest rates through time as a function of size of adjustment penalties, the latter being calculated as percentage of the maximum expected NPV in the basic case with no policy adjustment costs.", fig.height=3, fig.width=8----
Fig4a <- ggplot(harvest_trends, 
                aes(fraction, Ev, ymin=Ev-SDv, ymax=Ev+SDv, col=penalty)) + 
  geom_ribbon(aes(fill=penalty, col=NA), lwd=0, alpha=.05) + 
  geom_line() + 
  xlab("Fraction of NPV lost to costs") +
  scale_x_continuous(limits=c(0, 0.4)) + theme_tufte()


Fig4b <- ggplot(harvest_trends, 
                aes(fraction, Ea, ymin=Ea-SDa, ymax=Ea+SDa, col=penalty)) + 
  geom_ribbon(aes(fill=penalty, col=NA), lwd=0, alpha=.05) + 
  geom_line() + xlab("Fraction of NPV lost to costs")+
  scale_x_continuous(limits=c(0, 0.4)) + theme_tufte()


multiplot(Fig4a, Fig4b, cols=2)


## ----Figure_5, fig.cap = "(a) Variance of optimal stock sizes through time as a function of size of adjustment penalties, the latter being calculated as percentage of the maximum expected NPV in the basic case with no policy adjustment costs. (b) Autocorrelation of optimal stock sizes through time as a function of size of adjustment penalties, the latter being calculated as percentage of the maximum expected NPV in the basic case with no policy adjustment costs.", fig.height=3, fig.width=8----
Fig5a <- ggplot(fishstock_trends, 
                aes(fraction, Ev, ymin=Ev-SDv, ymax=Ev+SDv, col=penalty)) + 
  geom_ribbon(aes(fill=penalty, col=NA), lwd=0, alpha=.05) + 
  geom_line() + xlab("Fraction of NPV lost to costs")+
  scale_x_continuous(limits=c(0, 0.4))  + theme_tufte()


Fig5b <- ggplot(fishstock_trends, 
                aes(fraction, Ea, ymin=Ea-SDa, ymax=Ea+SDa, col=penalty)) + 
  geom_ribbon(aes(fill=penalty, col=NA), lwd=0, alpha=.05) + 
  geom_line() + xlab("Fraction of NPV lost to costs") +
  scale_x_continuous(limits=c(0, 0.4))  + theme_tufte()


multiplot(Fig5a, Fig5b, cols=2)



## ----Figure_1S, dependson="tidy", fig.cap="Multiple realizations of optimal harvesting strategy and resulting fish stock sizes under the different functional forms of adjustment costs. Note that while fish stock dynamics appear similar in penalty and non-penalty scenarios, the harvest dyanmics show clear patterns depending on the type of penalty.", fig.height=3, fig.width=8----

fig1s_df <- as.data.frame(subset(dt, replicate %in% paste0('rep_', 1:5)))
fig1s_df <- fig1s_df[c("time", "fishstock", "alternate", "harvest", "harvest_alt", "penalty_fn", "replicate")]
fig1s_df <- melt(fig1s_df, id = c("time", "penalty_fn", "replicate"))
fig1s_df <- data.frame(fig1s_df, baseline = fig1s_df$variable)
variable_map <- c(fishstock = "fish_stock", alternate = "fish_stock", harvest = "harvest", harvest_alt = "harvest")
baseline_map <- c(fishstock = "penalty", alternate = "no_penalty", harvest = "penalty", harvest_alt = "no_penalty")
fig1s_df$variable <- variable_map[fig1s_df$variable]
fig1s_df$baseline <- baseline_map[fig1s_df$baseline]


ggplot(fig1s_df, aes(time, value, col=baseline, lty=replicate)) +
  geom_line() +
  facet_grid(variable~penalty_fn) + 
  labs(x="time", y="stock size", title = "Example Stock & Harvest Dynamics")  + 
  theme_tufte()




