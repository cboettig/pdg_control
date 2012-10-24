#' SDP under multiple uncertainty
#'
#' Computes the SDP solution under the case of growth noise, 
#' implementation errors in harvesting, and meaurement errors in 
#' the stock assessment.  
#'
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma is the shape parameters for noise distribution (sigma_g, sigma_m, sigma_i) (default is no noise)
#' @param pdfn is the probability density function (same functional form is used for growth, measure, implement).  (Default is uniform)
#' @param profit is the profit function (defaults to the realized harvest)
#' @return The D matrix giving the optimal harvest for each possible state in each timestep.  
#' Hence the optimal policy at time t is given by the policy function D[,t].  
#' Values of D[,t] correspond to the index of h_grid.  Indices of of D[,t] correspond to states in y_grid.  
#' @export
SDP_multiple_uncertainty <- function(f, p, x_grid, h_grid, Tmax = 20,
                                     sigmas =c(sigma_g=0, sigma_m=0, sigma_i=0), 
                                     pdfn = unif_pdfn, profit = function(x,h) min(x, h)){
  
    ## Initialize stuff
    sigma_g <- sigmas['sigma_g']      
    sigma_m <- sigmas['sigma_m']
    sigma_i <- sigmas['sigma_i']
    D <- matrix(NA, nrow=length(x_grid), ncol=Tmax)  
        
    ## Compute the transition matrices 
    P <- outer(x_grid, h_grid, profit)
    F <- outer(x_grid, f(x_grid, 0, p), pdfn, sigma_g)
    M <- outer(x_grid, x_grid, pdfn, sigma_m)
    I <- outer(h_grid, h_grid, pdfn, sigma_i)
    
    # If no uncertainty in measure or implement
    V <- P
    for(t in 1:Tmax){
      D[,(Tmax-t+1)] <- apply(V, 1, which.max)
      V <- P + delta * F %*% V %*% D[, (Tmax-t+1)]  
    }
    
    # Uncertainty 
  #  Ep <- M %*% P %*% I
  #  U <- t(M) %*% F %*% M 
        
  #  V <- Ep
  #  for(t in 1:Tmax){
  #    D[,(Tmax-t+1)] <- apply(V, 1, which.max)
  #    V <- Ep + delta * U %*% V %*% I %*% D[, (Tmax-t+1)]  
  #  }
    
    D
}


unif_pdfn <- function(P, mu, s){
  if(mu==0){
    if(P == 0){
      1
    } else {
      0
  } else if(sigma > 0){
    dunif(P, mu*(1-s), mu*(1+s))
  } else { #sigma == 0
    bin_y <-  x_grid[which.min(abs(x_grid - P))]
    bin_fx <- x_grid[which.min(abs(x_grid - mu))]
    as.integer(bin_fx == bin_y)
  }
}
