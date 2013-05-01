#' Computes the transition matrix under the case of parametric uncertainty
#'
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param pardist a matrix with columns as variables and rows as the 
#' Monte Carlo samples from the posterior
#' @export
#' #FIXME: pardist should be in "long" format of par_posterior, allow for variable parameter numbers
parameter_uncertainty_SDP <- function(f, p, x_grid, h_grid, pardist, n_mc = 100){
  lapply(h_grid, function(h){ # for each h
    # Set up monte carlo sampling 
    d <- dim(pardist)
    indices <- round(runif(n_mc,1, d[1]))
    # Compute the transition matrix
    F_true <- 
      sapply(x_grid, function(x_t){           # For each x_t  
        bypar <- sapply(indices, function(i){ # For each parameter value
          p <- unname(pardist[i,c(2,1,3)])    # Parameters for the mean, at current sample 
          mu <- f(x_t,h,p)                    # expected x_{t+1} 
          est_sigma_g <- pardist[i,4]         # Variance parameter         
          # Calculate the probability density
          if(snap_to_grid(mu,x_grid) < x_grid[2]){ # handle the degenerate case 
            out <- numeric(length(x_grid))
            out[1] <- 1
            out
          } else { # Calculate the transition density based on log-normal
            # FIXME handle normally distributed measurement error additionally 
            out <- dlnorm(x_grid/mu, 0, est_sigma_g) 
          }
        })  
        # Average over parameter distributions and normalize
        ave_over_pars <- apply(bypar, 1, sum) 
        ave_over_pars / sum(ave_over_pars)
      })
    F_true <- t(F_true)
  })
}


snap_to_grid <- function(x, grid) sapply(x, function(x) grid[which.min(abs(grid - x))])   
