#' Determine the SDP matrix for uniform noise under multiple uncertainty
#'
#' Computes the transition matrix under the case of growth noise, 
#' implementation errors in harvesting, and meaurement errors in 
#' the stock assessment.  Assumes noise sources are distributed 
#' by uniform distributions to obtain analytic transition probability densities.
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
#' @param f_int is the function given by the analytic solution,  
#' $F(x,q) = \frac{1}{4\sigma_i\sigma_m} \int_{q-\sigma_i}^{q+\sigma_i} dh \int_{x-\sigma_m}^{x+\sigma_m} f(y, h)$
#' For logistic map $f(y,h) = 2 (y-h) \left(1-\frac{y-h}{2K}\right)$, this is: 
#' 
#' (where we've used m = $\sigma_m$, n = $\sigma_i$ for simplicity of Roman notation)
#' @export
SDP_uniform <- function(f, p, x_grid, h_grid, sigma_g,
                        pdfn=function(P, s) dunif(P, 1-s, 1+s),
                        sigma_m, sigma_i, f_int){
  gridsize <- length(x_grid)
  SDP_Mat <- lapply(h_grid, function(h){
    SDP_matrix <- matrix(0, nrow=gridsize, ncol=gridsize)
    # Cycle over x values
    for(i in 1:gridsize){ 
      ## Calculate the expected transition  
      x1 <- x_grid[i]
      x2_expected <- f_int(x1, h, sigma_m, sigma_i, p)
      ## If expected 0, go to 0 with probabilty 1
      if( x2_expected == 0) 
        SDP_matrix[i,1] <- 1  
      else {
        # relative probability of a transition to that state
        ProportionalChance <- x_grid / x2_expected
        Prob <- pdfn(ProportionalChance, sigma_g)

        # Store normalized probabilities in row
        SDP_matrix[i,] <- Prob/sum(Prob)
      }
    }
    SDP_matrix
  })
  SDP_Mat
}


#' Integrate against uniform distributions of measurement and implementation uncertainty
#' 
#' @export
int_f <- function(f, x, q, sigma_m, sigma_i, pars){
  K <- pars[2]
  sigma_m <- K*sigma_m
  sigma_i <- K*sigma_i # scale noise into units of K
  
  if(sigma_m > 0 && sigma_i > 0){
    g <- function(X) f(X[1], X[2], pars)
    lower <- c(max(x - sigma_m, 0), max(q - sigma_i, 0))
    upper <- c(x + sigma_m, q + sigma_i)
    A <- adaptIntegrate(g, lower, upper)
    out <- A$integral/((q+sigma_i-max(q-sigma_i, 0))*(x+sigma_m-max(x-sigma_m, 0)))
  } else if(sigma_m == 0 && sigma_i > 0){ 
    g <- function(h) f(x, h, pars)
    lower <- max(q - sigma_i, 0)
    upper <- q + sigma_i
    A <- adaptIntegrate(g, lower, upper)
    out <- A$integral/(q+sigma_i-max(q-sigma_i, 0))
  } else if(sigma_i == 0 && sigma_m > 0){ 
    g <- function(y) f(y, q, pars)
    lower <- max(x - sigma_m, 0)
    upper <- x + sigma_m
    A <- adaptIntegrate(g, lower, upper)
    out <- A$integral/(x+sigma_m-max(x-sigma_m, 0))
  } else if (m == 0 && n == 0){
    out <- f(x,q,pars)
  } else {
    stop("distribution widths cannot be negative")
  }
  out
}

#' Analytic solution to int_f for logistic model
#' @export
F_integral <- function(x,q, m, n, pars){
  K <- pars[2]
  m <- K*m
  n <- K*n # scale noise by K
  if(m > 0 && n > 0){
  out <- ((q+n-max(0,q-n))*(x+m-max(0,x-m))*(6*x*K-6*q*K-6*n*K+6*m*K+6*max(0,x-m)*K-6*
    max(0,q-n)*K-2*x^2+3*q*x+3*n*x-4*m*x-2*max(0,x-m)*x+3*max(0,q-n)*x-2*q^2-4*n*q+3*m*q+3*
    max(0,x-m)*q-2*max(0,q-n)*q-2*n^2+3*m*n+3*max(0,x-m)*n-2*max(0,q-n)*n-2*m^2-2*max(0,x-m)*m+3*
    max(0,q-n)*m-2*max(0,x-m)^2+3*max(0,q-n)*max(0,x-m)-2*max(0,q-n)^2))/(6*K)/((q+n-max(q-n, 0))*(x+m-max(x-m, 0)))
  } else if(m == 0 && n > 0) {
    y <- x
    out <- (((6*q+6*n-6*max(0,q-n))*y-3*q^2-6*n*q-3*n^2+3*max(0,q-n)^2)*K+(-3*q-3*n+3*max(0,q-n))*
      y^2+(3*q^2+6*n*q+3*n^2-3*max(0,q-n)^2)*y-q^3-3*n*q^2-3*n^2*q-n^3+max(0,q-n)^3)/(3*K*(q+n-max(q-n, 0)))
  } else if(n == 0 && m > 0){
    h <- q
    out <- ((3*x^2+(6*m-6*h)*x+3*m^2-6*h*m+6*max(0,x-m)*h-3*max(0,x-m)^2)*K-x^3+(3*h-3*m)*x^2+
      (-3*m^2+6*h*m-3*h^2)*x-m^3+3*h*m^2-3*h^2*m+3*max(0,x-m)*h^2-3*max(0,x-m)^2*h+max(0,x-m)^3)/(3*K*(x+m-max(x-m, 0)))
  } else if (m == 0 && n == 0){
    S <- max(x - q, 0)
    out <- S * (1 - S/K) + S
  } else {
    stop("distribution widths cannot be negative")
  }
  max(out,0)
}
