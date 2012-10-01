#' Determine the SDP matrix for uniform noise under multiple uncertainty
#'
#' Computes the transition matrix under the case of growth noise, 
#' implementation errors in harvesting, and meaurement errors in 
#' the stock assessment.  
#'
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma_g is the shape parameter (width) of the multiplicitive growth noise 
#' @param pdfn is the shape of the growth noise, which need not be uniform (is by default)
#' @param sigma_m is the half-width of the uniform measurement error (assumes uniform distribution)
#' @param sigma_i is the half-width of the implementation noise (always assumes uniform distribution)
#' @export
SDP_multiple_uncertainty <- function(f, p, x_grid, h_grid, sigma_g,
                        pdfn=function(P, mu, s){
                          if(mu==0)
                            if(P == 0)
                              1
                            else 
                              0
                          else 
                          dunif(P,mu*(1-s),mu*(1+s))
                          },
                        sigma_m, sigma_i){
  
    gridsize <- length(x_grid)
    
    transition <- function(x, y, sigma, g, pdfn){ # function for an entry of the matrix
      if(sigma > 0 )
        P <- pdfn(y, g(x), sigma)
      else {
        bin_y <-  x_grid[which.min(abs(x_grid - y))]
        bin_fx <- x_grid[which.min(abs(x_grid - g(x)))]
        P <- as.integer(bin_fx == bin_y)
      }
      P
    }
    
    m_grid <- expand.grid(x=x_grid, y=x_grid)
    F <- matrix(mapply(transition, m_grid$x, m_grid$y, 
                       MoreArgs = list(sigma=sigma_g, g=function(x) f(x,0,p), pdfn=pdfn)), nrow = gridsize)
    M <- matrix(mapply(transition, m_grid$x, m_grid$y,
                       MoreArgs = list(sigma=sigma_m, g=function(x) x, pdfn=pdfn)), nrow = gridsize)
    for(i in 1:gridsize){ # normalize
      F[i,] = F[i,]/sum(F[i,])
      M[i,] = M[i,]/sum(M[i,])
    }
    
  # Cycle over action space (harvest level)
  SDP_Mat <- lapply(h_grid, function(h){
    I <- matrix(mapply(transition, m_grid$x, m_grid$y, 
                       MoreArgs = list(sigma=sigma_i, g=function(x) max(x-h, 0), pdfn=pdfn)), nrow = gridsize)
    for(i in 1:gridsize)
      I[i,] = I[i,]/sum(I[i,])
    out <- F %*% M %*% I
    for(i in 1:gridsize)
      out[i,] = out[i,]/sum(out[i,])
    
    
    ## Debug Testing!
    F <- matrix(mapply(transition, m_grid$x, m_grid$y, 
                       MoreArgs = list(sigma=sigma_g, g=function(x) f(x,h,p), pdfn=pdfn)), nrow = gridsize)
    for(i in 1:gridsize) # normalize
      F[i,] = F[i,]/sum(F[i,])
    out <- F 
    
    
    out
  })
  SDP_Mat
}



