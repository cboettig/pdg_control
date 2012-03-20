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
<!--begin.rcode
p_grid = seq(0,.9,length=6)
x_grid = seq(1,10,length=7)
end.rcode-->


Define the two possible models:
<!--begin.rcode
sigma_g = 0.2

f1 = function(x_t1, x_t0){
  a = 1.5
  b = 0.05
  mu = a * x_t0 / (1 - b * x_t0)
  (mu <= 0) * (x_t1 == 0) +
  (mu > 0) * dlnorm(x_t1, log(mu) - sigma_g ^ 2 / 2, sigma_g)
}

f2 = function(x_t1, x_t0){
  A = 1
  B = 0.05
  mu = A * x_t0 / (1 - B * x_t0)
  (mu <= 0) * (x_t1 == 0) +
  (mu > 0) * dlnorm(x_t1, log(mu) - sigma_g ^ 2 / 2, sigma_g)
}
end.rcode-->


Define the transition probability function for going from any state `x_t0` and belief (i.e. probability that model 1 is true) `p_t0` to any other state/belief `x_t1`, `p_t1`.  Note that beliefs are updated by the simple Bayesian learning rule.
<!--begin.rcode
f = function(x_t0, p_t0, x_t1, p_t1){
  bw = (p_grid[2] - p_grid[1] ) / 2 # bin-width
  y1 = f1(x_t1, x_t0)
  y2 = f2(x_t1, x_t0)
  if(y1 == 0 && y2 == 0){
    P1 = p_t0
  } else {
    P1 = (p_t0 + bw) * y1
    P1 = P1 / (P1 + (1 - (p_t0+bw)) * f2(x_t1, x_t0))
    P1 = p_grid[p_grid < P1]
    P1 = P1[length(P1)]
 }
   y1 * ( p_t1 == P1)
}
end.rcode-->

A simple way to use this function to generate the matrix of all possible transitions (with thanks to [some SO folks](http://stackoverflow.com/questions/9652079/elegant-way-to-loop-over-a-function-for-a-transition-matrix-in-2-dimensions-in-r/9652497#9652497))
<!--begin.rcode
model_uncertainty <- function(x_grid, p_grid){
  bw = (p_grid[2] - p_grid[1] ) / 2 # bin-width
  d = expand.grid(x_t0 = x_grid, p_t0 = p_grid, x_t1 = x_grid, p_t1 = p_grid)
  M = matrix(mapply(f, d$x_t0, d$p_t0, d$x_t1, d$p_t1), nrow = length(p_grid) * length(x_grid) )
  for(i in 1:dim(M)[1])
    M[i,] = M[i,]/sum(M[i,])
  M
}
end.rcode-->
A modified version of finding the dynamic programming solution.  Not sure I've gotten this correct yet. 
<!--begin.rcode
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
      harvest_matrices(M, h_grid[i], x_grid, p_grid) %*% V + 
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
<!--begin.rcode
require(pdgControl)
M <- model_uncertainty(x_grid, p_grid)
h_grid <- seq(0,5, length=5)
T <- 5
xT <- 0
profit <- profit_harvest() 
delta <- 0.05
reward <- 0
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
end.rcode-->


<!--begin.rcode
bevholt <- function(x,h, p) p[1] * (x-h) / (1 - p[2] * (x-h))
sdp <- determine_SDP_matrix(bevholt, c(1.5, 0.05), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static$D
end.rcode-->

<!--begin.rcode
#extract_policy <- function(D, p_i, nx, np) D[p_i+(1:nx - 1) * np,]
extract_policy <- function(D, p_i, nx, np) D[(p_i-1)*nx + 1:nx,]
end.rcode-->

<!--begin.rcode
extract_policy(active$D, length(p_grid), length(x_grid), length(p_grid)) 
static$D
end.rcode-->



### Practice
Before writing the method to create the transition matrix for each value of the control variable (harvest) over the belief grid, I tested this simple example. 
<!--begin.rcode
## Assumes things are in blocks of shared x (i.e. vary over p, than vary over x)!
harvest_matrices = function(M, h, x_grid, p_grid){
  nx <- length(x_grid)
  np <- length(p_grid)
  harvest_all <- M[1:np,] ## assumes x_grid[1] = 0
  M_ <- M
  if(sum(x_grid < h) > 0){
    s <- which(x_grid < h)
    s <- s[length(s)] - 1 
    for(i in 1:nx){
      j = 1+(i-1)*np
      k = 1+(i-s-1)*np
      if(x_grid[i] < h)
        M_[j:(j+np-1),] <- harvest_all
      else
        M_[j:(j+np-1),] <- M[k:(k+np-1),]
    }
  }
  M_
}
end.rcode-->





We can generate the matrix for each `h` value from the rows of the existing matrix, rather than editing all the functions above to depend on `h` as well and looping over all possible harvest values. Requires our h_grid equal our x_grid 
<!--begin.rcode
## Assumes things are in blocks of shared p (i.e. vary over x, than vary over p)!
harvest_matrices = function(M, h, x_grid, p_grid){
  nx <- length(x_grid)
  np <- length(p_grid)
  # x is the integer index along x_grid, not the floating pt value
  extract_fixed_x <- function(x) x + (1:np - 1) * np
  extract_fixed_p <- function(p) (p-1)*nx + 1:nx

  dead <- numeric(nx*np)
  dead[1+(1:np-1)*nx] <- 1
  kill <- sapply(which(x_grid < h), extract_fixed_x)
  M[kill,] <- dead
  
  

  M[harvested,] <- M[unharvested,]

  harvest_all <- M[1:np,] ## get all x_grid[1] = 0
  M_ <- M
  if(sum(x_grid < h) > 0){
    s <- which(x_grid < h)
    s <- s[length(s)] - 1 
    for(i in 1:nx){
      j = 1+(i-1)*np
      k = 1+(i-s-1)*np
      if(x_grid[i] < h)
        M_[j:(j+np-1),] <- harvest_all
      else
        M_[j:(j+np-1),] <- M[k:(k+np-1),]
    }
  }
  M_
}
end.rcode-->



<!--begin.rcode
matrix_given_harvest = function(M, h, x_grid){
  M_ <- M
  zeros <- x_grid < h
  if(sum(zeros) != 0){
    harvest_all <- numeric(length(x_grid))
    harvest_all[1] <- 1
    M_[zeros,] <- matrix(rep(harvest_all, sum(zeros)), nrow=sum(zeros), byrow=T)
    M_[!zeros,] <- M[1:sum(!zeros),]
  }
  M_
}
end.rcode-->
