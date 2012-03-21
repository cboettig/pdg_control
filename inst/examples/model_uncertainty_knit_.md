<!--roptions dev="png", fig.width=7, fig.height=5, fig.path='ex-out-', tidy=FALSE, warning=FALSE, comment=NA, message=FALSE, cache=FALSE-->

<!--begin.rcode echo=FALSE 
render_gfm()
opts_knit$set(upload = TRUE)
## use flickr to upload with these options
require(socialR)
options(flickrOptions=list(
  description="https://github.com/cboettig/pdg_control/blob/master/inst/examples/",
  tags="stochpop, pdg_control"))
opts_knit$set(upload.fun = flickr.url)
end.rcode-->


## Model Uncertainty
Notes from my first attempt at coding the active adaptive management solution to the simple model-uncertainty problem. 


Set a moderate example grid
<!--begin.rcode grids
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(1,10,length=10) 
sigma_g = 0.2
end.rcode-->

Define the transition densities for two different models:
<!--begin.rcode f1
f1 = function(x_t1, x_t0){
  a = 1.5
  b = 0.05
  mu = a * x_t0 / (1 - b * x_t0)
  (mu <= 0) * (x_t1 == 0) +
  (mu > 0) * dlnorm(x_t1, log(mu), sigma_g)
}
end.rcode-->

<!--begin.rcode f2
f2 = function(x_t1, x_t0){
  A = 1.1
  B = 0.05
  mu = A * x_t0 / (1 - B * x_t0)
  (mu <= 0) * (x_t1 == 0) +
  (mu > 0) * dlnorm(x_t1, log(mu), sigma_g)
}
end.rcode-->


Define the transition probability function for going from any state `x_t0` and belief (i.e. probability that model 1 is true) `p_t0` to any other state/belief `x_t1`, `p_t1`.  Note that beliefs are updated by the simple Bayesian learning rule.
<!--begin.rcode learning
f = function(x_t0, p_t0, x_t1, p_t1){
  y1 = f1(x_t1, x_t0)
  y2 = f2(x_t1, x_t0)
  P1 = p_t0 * y1
  P1 = P1 / (P1 + (1 - p_t0) * y2)
  if(is.na(P1))
    P1 = p_t0
  else{
    i <- 1
    np <- length(p_grid)
    while(p_grid[i] < P1 & i < np)
      i <- i+1
    P1 <- p_grid[i]  
 }
 y1 * ( p_t1 == P1)
}

end.rcode-->

<!--begin.rcode matrix
# f = function(a,b,c,d) sprintf("%.1f, %.1f, %.1f, %.1f", a,b,c,d)

model_uncertainty <- function(x_grid, p_grid, h_grid){
  Matrices <- lapply(h_grid, function(h){
    x_minus_h <- (x_grid-h) * as.integer( (x_grid-h)>0 )

    d = expand.grid(x_t0 = x_minus_h, p_t0 = p_grid, x_t1 = x_grid, p_t1 = p_grid)
    M = matrix(mapply(f, d$x_t0, d$p_t0, d$x_t1, d$p_t1), nrow = length(p_grid) * length(x_grid) )

## old way
#    M <- matrix(
#    sapply(p_grid, function(p_t1)
#      sapply(x_grid, function(x_t1)
#        sapply(p_grid, function(p_t0)
#          sapply(x_minus_h, function(x_t0)
#            f(x_t0, p_t0, x_t1, p_t1) )))),
#    nrow = length(p_grid) * length(x_grid))

    
    for(i in 1:dim(M)[1]) # normalize
      M[i,] = M[i,]/sum(M[i,])
    M })
  Matrices
}
end.rcode-->


A modified version of finding the dynamic programming solution.  Not sure I've gotten this correct yet. 
<!--begin.rcode define_dp
dp_optim <- function(M, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0, p_grid){
  gridsize <- length(x_grid) * length(p_grid)
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  profit.grid <- function(x_grid, h_i)
    expand.grid(profit(x_grid, h_i), p_grid)[[1]]

  # give a fixed reward for having value larger than xT at the end. 
  V[1+(x_grid-1)*length(p_grid) >= xT] <- reward # a "scrap value" for x(T) >= xT

  # loop through time  
  for(time in 1:OptTime){ 
    # try all potential havest rates
    V1 <- sapply(1:HL, function(i){
      # Transition matrix times V gives dist in next time
      M[[i]] %*% V + 
      # then (add) harvested amount times discount
       profit.grid(x_grid, h_grid[i]) * (1 - delta) 
    })

    # find havest, h that gives the maximum value
    out <- sapply(1:gridsize, function(j){
      value <- max(V1[j,], na.rm = T) # each col is a diff h, max over these
      index <- which.max(V1[j,])  # store index so we can recover h's 
      c(value, index) # returns both profit value & index of optimal h.  
    })
    # Sets V[t+1] = max_h V[t] at each possible state value, x
    V <- out[1,]                        # The new value-to-go
    D[,OptTime-time+1] <- out[2,]       # The index positions
  }
  # Format the output 
  list(D=D, V=V)
}
end.rcode-->

Sticking the pieces together,
<!--begin.rcode pars
require(pdgControl)
h_grid <- x_grid-1 
T <- 5
xT <- 0
profit <- profit_harvest(price=10, c0=30) 
delta <- 0.05
reward <- 0
end.rcode-->

Active Adaptive Mangement solution
<!--begin.rcode active
M <- model_uncertainty(x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
end.rcode-->

Static solution
<!--begin.rcode static
bevholt <- function(x,h, p) p[1] * (x-h) / (1 - p[2] * (x-h))
sdp <- determine_SDP_matrix(bevholt, c(1.5, 0.05), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static$D
end.rcode-->

Utils
<!--begin.rcode utils
indices_fixed_x <- function(x, np, nx) (1:np-1)*nx + x
indices_fixed_p <- function(p, nx) (p-1)*nx + 1:nx
extract_policy <- function(D, p_i, nx, np) D[(p_i-1)*nx + 1:nx,]
end.rcode-->


Confirm that the policy with high probability on model 1 matches the static solution for model 1:
<!--begin.rcode
extract_policy(active$D, length(p_grid), length(x_grid), length(p_grid)) 
static$D
end.rcode-->









Another sanity check
<!--begin.rcode
np <- length(p_grid)
nx <- length(x_grid)
sdp2 <-
lapply(1:length(M), function(i){
  collapse_to <- sapply(1:nx, function(i) rowSums(M[[i]][,indices_fixed_x(i,np,nx)]))
  collapse_to[(nx*np-nx+1):(np*nx),]
})
find_dp_optim(sdp2, x_grid, h_grid, T, xT=0, profit, delta, reward)$D
end.rcode-->


