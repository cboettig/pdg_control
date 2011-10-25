# file: training_prob2_collocation.R
# author: Carl Boettiger http://carlboettiger.info
# date: 2011-11-24

# state equation
alpha <- 1
K <- 1
C <- 0.2
gamma <- -2
rho <- 0.05

f <- function(x, h, t) 
  x * alpha * ((K - x) / K) * ((x - C) / K) - h*x
df <- function(x, h, t) 
  - alpha * ( C* K - 2 * x * K - 2 * x * C + 3 * x ^ 2) / K ^ 2 - h

#f <- function(x, h, t) 
#  x * alpha * (1 - x / K) - h * x
#df <- function(x, h, t) 
#  alpha * (1 - x / K) - alpha * x / K 

require(deSolve)
fish <- function(t,y,p){
  yd <- f(y[1], p[1], t)
  list(c(yd), NULL)
}
t <- seq(0, 10, by=0.1)
p <- c(0.0065)
out <- lsoda(0.9, t, fish, p, rtol=1e-4)
plot(out[,1], out[,2], type="l")


require(bvpSolve)
pars <- c(gamma, K, C, rho)

#' @param x variable.  y' = f( y(x) )
#' @param y a vector of the system state y[1],y[2] = (x,h) 
#' y1' = f(y1, y2)
#' y2' =  \gamma^{-1} (\rho - f'(x) )
#' x \in (0, 1), y(x_0) = (1000,NA), y(x_T) = (1000, NA)
fun <- function(x,y,pars){
  gamma <- pars[1]
  K <- pars[2]
  C <- pars[3]
  rho <- pars[4]
  dy1 <- f(y[1], y[2], x) 
  dy2 <- (1 / gamma) * (rho - df(y[1], y[2], x))
  list(c(dy1, dy2))
}

x <- seq(0, 10, by=1)
sol1 <- bvpshoot(yini = c(1,NA), yend = c(0.35, NA),
                 x = x, func = fun, guess = 0, parms=pars)

sol3 <- bvpcol(yini = c(1, NA), yend = c(.35, NA),
               x = x, func = fun, parms=pars, 
               xguess = sol1[,1], yguess=t(sol1[,-1]))

plot(sol3[,1], sol3[,2], type="l", lwd=3, ylim=c(0,1), xlab="time", ylab="fraction of K")
lines(sol3[,1], sol3[,3], col="red", lwd=3)
legend("topright", c("fish stock", "harvest effort"), lty=1, col=c("black", "red"))

diagnostics(sol3)

png("collocation.png")
  plot(sol3[,1], sol3[,2], type="l", lwd=3, ylim=c(0,1), xlab="time", ylab="fraction of K")
  lines(sol3[,1], sol3[,3], col="red", lwd=3)
  legend("topright", c("fish stock", "harvest effort"), lty=1, col=c("black", "red"))
dev.off()


# sol2 <- bvptwp(yini = c(1, NA), yend = c(.35, NA),
#               x = x, func = fun, parms=pars)



