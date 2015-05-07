parallel = FALSE
ncpu <- 4

library("pdgControl")
library("tidyr")
library("dplyr")
library("ggplot2")
library("reshape2")
library("data.table")
library("snowfall")

## ----profit_model--------------------------------------------------------
price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)
profit_quad <- function(c2) profit_harvest(price = price, c0 = c0, c1 = c2)


## ----c2_grid-------------------------------------------------------------
c2 <- exp(seq(0, log(41), length.out = 40))-1
c2 <- seq(0, 40, length.out=100)
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
none <- function(h, h_prev)  0
i <- which(x_grid > K)[1]

sfInit(parallel = parallel, cpus = ncpu)
sfLibrary(pdgControl)
sfExportAll()
policies <- sfLapply(c2, function(c2){
  policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                         profit_quad(c2), delta, reward, penalty = none)
})


fee <- sapply(policies, function(policy) max(policy$V[i,]))
npv0 <- fee[1]
normalized_fee <- (npv0-fee)/npv0

apples_index <- which.min(abs(normalized_fee - reduction))
apples <- c2[apples_index]

quad_policy <- policies[[apples_index]]$D

## ----simulate_policy, dependson=c("policynames", "apples")---------------
reps <- 1:100
names(reps) = paste("rep", 1:length(reps), sep="_") # treat as a factor
seeds <- 1:100
sims <- lapply(reps, function(x) 
                simulate_optim(f, pars, x_grid, h_grid, x0, 
                               quad_policy, z_g, z_m, z_i, 
                               opt$D, profit=profit_quad(c2 = apples), 
                               penalty=none, seed=seeds[x]))


## ----tidy, dependson="simulate_policy"-----------------------------------
#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]]))  
dt_quad <- data.table(dat)
setnames(dt_quad, "L1", "replicate") # names are nice



fraction_no_shift <- function(x){
  noshift <- x[2:length(x)] - x[1:length(x)-1] == 0
  strictly_positive <- (x > 0.01 )[-1] 
  who <- noshift & strictly_positive 
  sum(who)/length(who)
}

