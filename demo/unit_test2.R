require(PDGcontrol)

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
sigma_g <- 0.2      # Noise in population growth
gridsize <- 10   # gridsize (discretized population)
sigma_m <- .2  # 
sigma_i <- .2 # 
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
h <- x_grid[3]
y <- x_grid[4]
require(cubature)

  expected <- f(y,h,p)
      if(expected==0){
        Prob <- numeric(gridsize)
        Prob[1] <- 1
      } else {
        # dividing x by the expected value is same as scaling distribution to mean 1
        pdf_zg <- function(x, expected) dlnorm(x/expected, 0, sigma_g)
        pdf_zm <- function(x) dlnorm(x, 0, sigma_m)
        pdf_zi <- function(x,q) dlnorm(x, log(q), sigma_i)
        Prob <- sapply(x_grid, function(y){
          F <- function(x) 
            pdf_zg(y, f(x[1], x[2], p)) * pdf_zm(x[1]) * pdf_zi(x[2], h)
          int <- adaptIntegrate(F, c(0, 0), c(10*K, 10*K))
          int$integral
        })   
      }
      Prob/sum(Prob)
