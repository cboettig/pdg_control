require(pdgControl)

sigma_g <- 0.8    # Noise in population growth
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
f <- BevHolt                # Select the state equation
pars <- c(2, 4)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity 


x <- numeric(1000)
x[1] <- .01
for(t in 1:(length(x)-1)){
  x[t+1] = z_g()*f(x,0,pars)
}

likelihood_estimator <- function(f, x, guess_p = NULL){
  loglikfn <- function(pars){
    probs <- sapply(1:(length(x)-1), function(i){
      mu <- f(x[i],0,pars)
      sigma <-  
      dlnorm(x[i+1], meanlog=log(mu+sigma_g^2/2), sdlog=sigma_g, log=TRUE)
    })
    -sum(probs)
  }
  o <- optim(guess_p, loglikfn)
}

out <- likelihood_estimator(f, x, pars)
out


X <- 5 * rlnorm(10000, 0, 1)
-sum(dlnorm(X, log(5), 1))



  loglik_noise <- function(pars){
    probs <- sapply(1:(length(x)-1), function(i){
      p <- pars[-1]    # estimate noise parameter too?
      sigma <- pars[1]
      dlnorm(x[i+1], logmean=log(f(x[i],0,p)), sdlog=log(sigma), log=TRUE)
    })
    -sum(probs)
  }

