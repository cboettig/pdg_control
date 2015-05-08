profit_quad <- function(c2) profit_harvest(price = price, c0 = c0, c1 = c2)

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

