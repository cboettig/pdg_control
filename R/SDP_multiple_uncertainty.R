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
SDP_multiple_uncertainty <- function(f, p, x_grid, h_grid, Tmax = 25,
                                     sigmas =c(sigma_g=0, sigma_m=0, sigma_i=0), 
                                     pdfn = pdfn, profit = function(x,h) pmin(x, h)){
  
    ## Initialize stuff
    sigma_g <- sigmas['sigma_g']      
    sigma_m <- sigmas['sigma_m']
    sigma_i <- sigmas['sigma_i']
    n_x <- length(x_grid)
    n_h <- length(h_grid)
    D <- matrix(NA, nrow=n_x, ncol=Tmax)  
            
    ## Compute the transition matrices 
    P <-  outer(x_grid, h_grid, profit) 
    M <- rownorm( outer(x_grid, x_grid, pdfn, sigma_m) )
    I <- rownorm( outer(h_grid, h_grid, pdfn, sigma_i) )
    
    # integrate over uncertainty in I and X 
    F <- lapply(1:n_h, function(q){  
      t(sapply(1:n_x, function(y){
        out <- numeric(n_x)
        mu <- sum(sapply(1:n_x, function(x)
          f(x_grid[x],h_grid,p) %*% 
            I[q,] *  # Implementation error
            M[x,y])) # Measurement error
        if(mu==0)
          out[1] <- 1
        else 
          out <- dlnorm(x_grid/mu, 0, sigma_g)
        out/sum(out)
      }))
    })
    
    Ep <- M %*% P %*% I
    V <- Ep
    for(t in 1:Tmax){
      D[,(Tmax-t+1)] <- apply(V, 1, which.max)
      v_t <- apply(V, 1, max) # vector of current values
      V <- sapply(1:n_h, function(j){ # updated value matrix
        Ep[,j] + (1-delta) * M %*% F[[j]] %*% v_t
      })
    }
    D
}
  
# row-normalize the probability distribution
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))

# binning routine (vectorized)
snap_to_grid <- function(x, grid) sapply(x, function(x) grid[which.min(abs(grid - x))])   

# Uniform pdf where width scales with mean (probably better way to do that)
# And handles degenerate/delta fn cases such as no width and zero mean
FUN <- function(P, mu, s){
  if(mu == 0){
    as.integer(P == 0)
  } else if(s > 0){
    if(mu > 0) dunif(P, mu * (1 - s), mu * (1 + s))
  }
  else { # delta spike
    P <- snap_to_grid(P, x_grid)
    mu <- snap_to_grid(mu, x_grid)
    as.integer(P == mu)
  }
}
pdfn <- Vectorize(FUN)


############## DEPRECATED SCRATCH ###########

dummy_attempts <- function(){
    #F <- rownorm( outer(x_grid, f(x_grid, 0, p), pdfn, sigma_g) ) # has fatal cases where there are values of f(x) for which no x can get within +/- f(x)*sigma_g of  
    ## This is  the F without uncertainty in current state (measurement, growth).  Alternately, use this F[[1]] (e.g. h=0) and the Q matrix below
    F <- lapply(h_grid, function(h){
      t(sapply(f(x_grid,h,p), function(y){
        out <- numeric(n_x)
        if(y==0)
          out[1] <- 1
        else 
         out <- dlnorm(x_grid/y, 0, sigma_g)
        out/sum(out)
      }))
    })
    F2 <- determine_SDP_matrix(f, p, x_grid, h_grid, sigma_g, pdfn=function(P, s) dlnorm(P, 0, s))  
    ## Without uncertainty in measure or growth, here's the SDP algorithm:
    V <- P
    for(t in 1:Tmax){
      D[,(Tmax-t+1)] <- apply(V, 1, which.max)
      v_t <- apply(V, 1, max) # vector of current values
      V <- sapply(1:n_h, function(j){ # updated value matrix
        P[,j] + (1-delta) * F[[j]] %*% v_t
      })
    }
    D   
    ## Q  is the map that multiplies F[[1]] to account for implmentation uncertainty.  
    ## This calculation is super slow.  In principle should work but is not matching up.  Probably need to transpose something to correct for the sapplys
    Q <- lapply(1:n_h, function(q)
      sapply(1:n_x, function(s)
        sapply(1:n_x, function(x)
          sum(sapply(1:n_h, function(h) 
            (x_grid[s] == pmax(x_grid[x]-h_grid[h],0)) * I[q,h]
    )))))
   
}

