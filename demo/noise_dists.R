#########################################################################
# A library of distribution functions for sources of stohcasticity      #
#########################################################################

# z_g is the stochasticity in the growth process, x_{t+1} = z_g f(x_t)
# z_m is the measurement error in the stock assessment, m_t = z_m x_t
# z_i is the implementation error in the harvest quota: h_t = z_i q_t

# Normal random vars -- an unusual choice given the negative domain support
#z_g <- function() rnorm(1,1, sigma_g)
#z_m <- function() rnorm(1,1, sigma_m)
#z_i <- function() rnorm(1,1, sigma_i)

# Log-normal distribution -- perhaps the most natural, at least for z_g
# mean is 1 = exp(mu + sigma^2/2), then
# log(1) - sigma^2/2 = mu
#z_g <- function() rlnorm(1,  log(1)-sigma_g^2/2, sigma_g) # mean 1
#z_m <- function() rlnorm(1,  log(1)-sigma_m^2/2, sigma_m) # mean 1
#z_i <- function() rlnorm(1,  log(1)-sigma_i^2/2, sigma_i) # mean 1

z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1

# Uniform distribution
#z_g <- function() runif(1, max(0,1-sigma_g), 1+sigma_g)
#z_m <- function() runif(1, max(0,1-sigma_m), 1+sigma_m)
#z_i <- function() runif(1, max(0,1-sigma_i), 1+sigma_i)


# No noise
#z_g <- function() 1
#z_m <- function() 1
#z_i <- function() 1


