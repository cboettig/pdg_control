
source("stochastic_dynamic_programming.R")

source("population_models.R")

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
sigma_g <- 0.2      # Noise in population growth
gridsize <- 100   # gridsize (discretized population)
sigma_m <- .0  # 
sigma_i <- .0 # 
interval <- 1

# Chose the state equation / population dynamics function
f <- BevHolt 
pars <- c(2,4)  
K <- (pars[1]-1)/pars[2] 
xT <- 0 
e_star <- 0

# Choose grids 
x_grid <- seq(0, 2*K, length=gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same res as stock


p <- pars
i <- 3
bw <- x_grid[2]-x_grid[1]
q <- x_grid[3]


require(cubature)
pdf_zg <- function(x, expected) dlnorm(x, log(expected)-sigma_g^2/2, sigma_g)
pdf_zm <- function(x) dlnorm(x, log(1)-sigma_m^2/2, sigma_m)
pdf_zi <- function(x,q) dlnorm(x, log(q)-sigma_i^2/2, sigma_i)
F <- function(x) pdf_zg(x[1], f(x[2],x[3],p))*pdf_zm(x[2])*pdf_zi(x[3], q)
adaptIntegrate(F, c(x_grid[i]-bw,0,0), c(x_grid[i]+bw, x_grid[gridsize], h_grid[gridsize]))
