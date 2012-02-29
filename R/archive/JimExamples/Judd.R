# file Judd.R
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-21
#  
# 
# Use the collocation method to solve boundary value problems
# This code solves the problem on page 390 in Judd
# which is a life-cycle consumption and savings model
#
#

### Definition of the problem, from page 353 (10.6.10) ###
# 
# Maximize the objective function:
# \[ \max_c \int_0^T e^{-\rho t} u(c) dt \]
# s.t. [\dot A = f(A) + w(t) -c(t) \]
# \[ A(0) = A(T) = 0 \]
# 

# with w(t) wage rate at time t
# A(t) assets at ime t, 
# and f(A) the return on invested assets

# From the statement of the problem, we have:
# the Hamiltonian: \( H = u(c) + \lambda ( f(A) + w(t) -c ) \)
# and the costate eqn: \( \dot \lambda =  \lambda (\rho  - f'(A) ) \)
# and bdry conditions.  

# assuming u concave in c, the maximum principle gives us:
# the first order condition that:
# 0 = \partial_c H = u'(c) - \lambda


# So we're left with a nice ODE system:
# \[ \dot A = f(A) + w(t) - C(\lambda) \]
# \[ \dot \lambda = \lambda (\rho - f'(A)) \]
# A(0) = A(T) = 0

# Let's assume functional forms for the utility:
# \[ u(c) = c^{1+\gamma}/(1+\gamma) \]
# and the asset return:
# \[ f(A) = rA \]

# Then our system becomes 
# \[ \dot A = rA + w(t) - c(t) \]
# \[ \dot \lambda = \lambda (\rho - r ) \]

# From the max principle of Hamiltonian, first order condition
# \[ 0 = \partial_c H \], so
# \[ u'(c) -\lambda = c^{\gamma} -\lambda \]
# \[ \lambda  = c^{\gamma} \]
# \[ \frac{d\lambda}{dt}  = \frac{d\lambda}{dc}\frac{d c}{dt} \] 
# \[ \frac{d\lambda}{dc} = \gamma c^{\gamma} \]
# which we sub into the costate equation to get
# \[ \frac{d\lambda}{dt} =   c^{\gamma} (\rho - r) \]
# \[ \gamma c^\gamma \frac{dc}{dt}  = c^{\gamma} (\rho - r) \]
# which simplfiies to: 
# \[ \frac{dc}{dt}  = \frac{c}{\gamma} (\rho - r) \]

# and our second equation is still:
# \[ \dot A = rA + w(t) - c(t) \]


# We can write these two equations as residuals to minimize:
# \[ R_1(t) = \dot c(t) - (\rho-r) c(t)/\gamma \]
# \[ R_2(t) = \dot A(t) - \left(r A(t) +  w(t)  - c(t) \right) \]

# let's further assume a functional form of w(t):
w <- function(t) sapply(t, 0.5 + t/10-4*(t/50)^2)



### Let's also define the Chebychev basis functions ####

#' Generate the Chebychev basis functions at the nodes
#' \[ T_n(x) = cos( n arccos x) \]
Chebybasis <- function(n,m,x){
    sapply(1:m, function(i){
      sapply(1:n, function(j) cos((j-1)*acos(x[i])))
    })
}
#' Generate the derivatives of the Cheby basis functions
Dchebybasis <- function(n, m, x){
  sapply(1:length(x), function(i){
    sapply(1:n, function(j) 
      (j-1)*sin((j-1)*acos(x[i]))/((1-x[i]^2)^.5) 
    )
  })
}
    


# Economic parameters
R <- 0.1      # a constant of the asset function f(A)
gam <- -2     # a constant of the utility function, u(c)
rho1 <- 0.05  # discount rate
W <- 1        # fixed wage (option)

# Boundary conditions
A0 <- 0 
AF <- 0

# Time span
t0 <- 0 
tf <- 50

# let's take a look at w(t) 
plot(1:tf, w(1:tf), type="l", main="Wage profile over time")

# Chebyshev nodes
n <- 11
m <- n
a <- t0
b <- tf



# what the polynomials look like:
plot(seq(-1,1,length=m), T[1,], ylim=c(-1,1), type="l")
for(i in 2:11) lines(seq(-1,1,length=m), T[i,], col=i, lwd=11-i)

## Why?
xn <-  sapply(1:m, function(i) cos(((m-i+0.5)/m)*pi))
## Why?  
Time <- sapply(1:m, function(i) (a+b)/2 + xn[i]*(b-a)/2)

# T_i(x) is the degree i Chebyshev polynomial:
# Generate T(x) functions at the nodes with the coefficients
T = Chebybasis(n,m,xn)
dT = (2 / (b-a) ) * Dchebybasis(n, m, xn)



c <- matrix(0, nrow=n, ncol=2)
c[1,] <- 1 # initial guess

options <- optimset



#' residuals function
#' @param c Control variable (consumption)
#' @param R constant of the asset function
#' @param gam constant of the utility funciton
#' @param rho1 discounting
#' @param W wage (if constant)
#' @param A0 boundary condition of assets at t0
#' @param AF bc of assets at end, tf
#' @param n discretization scale
#' @param a left Chebyschev node interval, t0 
#' @param b right interval, tf
#' @param xn Chebychev nodes
#' @param Time (?) from xn
#' @param T Chebychev polynomials
#' @param dT derivative of Chebychev polynomails
#' @return 
#' @details from Judd p390 (derivation above)
#' \[ R_1(t) = \dot c(t) - 0.025 c(t) \]
#' \[ R_2(t) = \dot A(t) - \left(0.1 A(t) +  w(t)  - c(t) \right) \]
resid1 <- function(c, R, gam, rho1, W, A0, AF, n, 
                   a, b, xn, Time, T, dT){
  co <- numeric(n)
  A <- numeric(n)

  co <- T %*% c[1,]
  A <- T %*% c[2,]
  cdot <- dT %*% c[1,]
  Adot <- dT %*% c[2,]

  R1 <- cdot - (1 / gam) * co * (rho1 - R)
  R2 <- Adot - (R * A + w(Time) -co)
### Fixed wage option:
# R2 <- Adt - R*A - W + co

  Res1 <- A0-c[,2] %*% (-1)^(1:n)
}



