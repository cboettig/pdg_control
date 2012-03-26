





## Model Uncertainty
Notes from my first attempt at coding the active adaptive management solution to the simple model-uncertainty problem. 


Set a moderate example grid


```r
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(1,10,length=10) 
sigma_g = 0.2
```




Define some utilities to handle the combined state-space/belief-space, `(x,p)`. 


```r
nx <- length(x_grid)
np <- length(p_grid)
indices_fixed_x <- function(x) (1:np-1)*nx + x
indices_fixed_p <- function(p) (p-1)*nx + 1:nx
extract_policy <- function(D, p_i, nx, np) D[(p_i-1)*nx + 1:nx,]
```





Define the transition densities for two different models:


```r
f1 = function(x_t1, x_t0){
  a = 1.5
  b = 0.05
  mu = a * x_t0 / (1 - b * x_t0)
  (mu <= 0) * (x_t1 == 0) +
  (mu > 0) * dlnorm(x_t1, log(mu), sigma_g)
}
```






```r
f2 = function(x_t1, x_t0){
   mu = 3 * x_t0 / (1 - 0.05 * x_t0)
  (mu <= 0) * (x_t1 == 0) +
  (mu > 0) * dlnorm(x_t1, log(mu), sigma_g)
}
```





Define the transition probability function for going from any state `x_t0` and belief (i.e. probability that model 1 is true) `p_t0` to any other state/belief `x_t1`, `p_t1`.  Note that beliefs are updated by the simple Bayesian learning rule.


```r
f = function(x_t0, p_t0, x_t1, p_t1){
  y1 = p_t0 * f1(x_t1, x_t0)
  y2 = (1-p_t0) * f2(x_t1, x_t0)
  P1 = y1 / (y1 + y2)
  if(is.na(P1) || x_t0 == 0)
    P1 = p_t0
  else{
    i = 1
    np = length(p_grid)
    while(p_grid[i] < P1 & i < np)
      i = i+1
    P1 = p_grid[i]  
 }
 (y1+y2) * ( p_t1 == P1)
}
```





Some unit tests of this behavior: Should you see a transition from 1 to 10, you should be almost sure it came from model 2, and hence move to the first bin where belief in model 1 is <!--inline.rcode p_grid[1]-->, even if you were 0.99 sure that model 1 was correct until then.


```r
sapply(p_grid, function(p) f(1,.99,10,p))
```



```
[1] 1.222e-10 0.000e+00 0.000e+00 0.000e+00 0.000e+00
```



Should you see a transition from 1 to 2, you should become almost sure model 1 is correct, even if it had only a 1% probability previously:


```r
sapply(p_grid, function(p) f(1,.01,2,p))
```



```
[1] 0.00000 0.07772 0.00000 0.00000 0.00000
```




But we want the amount you change your belief to depend on where you started.  These should be different:


```r
sapply(p_grid, function(p) f(1,.01,2,p))
```



```
[1] 0.00000 0.07772 0.00000 0.00000 0.00000
```



```r
sapply(p_grid, function(p) f(1,.99,2,p))
```



```
[1] 0.0000 0.0000 0.0000 0.0000 0.4918
```






A simple way to use this function to generate the matrix of all possible transitions (with thanks to [some SO folks](http://stackoverflow.com/questions/9652079/elegant-way-to-loop-over-a-function-for-a-transition-matrix-in-2-dimensions-in-r/9652497#9652497))



```r
model_uncertainty <- function(x_grid, p_grid, h_grid){
  lapply(h_grid, function(h){
    x_minus_h <- (x_grid-h) * as.integer( (x_grid-h)>0 )
    d = expand.grid(x_t0 = x_minus_h, p_t0 = p_grid, x_t1 = x_grid, p_t1 = p_grid)
    M = matrix(mapply(f, d$x_t0, d$p_t0, d$x_t1, d$p_t1), nrow = length(p_grid) * length(x_grid) )
    for(i in 1:dim(M)[1]) # normalize
      M[i,] = M[i,]/sum(M[i,])
    M
  })
}
```






A modified version of finding the dynamic programming solution.  Not sure I've gotten this correct yet. 


```r
dp_optim <- function(M, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0, p_grid){
  gridsize <- length(x_grid) * length(p_grid)
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  profit.grid <- function(x_grid, h_i)
    expand.grid(profit(x_grid, h_i), p_grid)[[1]]

  # give a fixed reward for having value larger than xT at the end. 
  V[sapply(x_grid[x_grid>xT], function(x) indices_fixed_x(x))] <- reward

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
```






Sticking the pieces together,


```r
require(pdgControl)
h_grid <- x_grid-1 
T <- 5
xT <- 0
profit <- profit_harvest(price=10, c0=30) 
delta <- 0.05
reward <- 0
```




Active Adaptive Mangement solution


```r
M <- model_uncertainty(x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
```





Let's make sure the matrix is working correctly.  Transitions from 1 to 2 should be going to the far right bins, representing model 1, while those from 1 to 10 should go to the far left, representing no faith in model 1. 


```r
M[[1]][indices_fixed_x(1), indices_fixed_x(2)]
```



```
     [,1]    [,2] [,3]   [,4]   [,5]
[1,]    0 0.07836    0 0.0000 0.0000
[2,]    0 0.00000    0 0.1999 0.0000
[3,]    0 0.00000    0 0.0000 0.3468
[4,]    0 0.00000    0 0.0000 0.5277
[5,]    0 0.00000    0 0.0000 0.7562
```



```r
M[[1]][indices_fixed_x(1), indices_fixed_x(10)]
```



```
          [,1] [,2] [,3] [,4] [,5]
[1,] 1.219e-08    0    0    0    0
[2,] 1.004e-08    0    0    0    0
[3,] 7.439e-09    0    0    0    0
[4,] 4.234e-09    0    0    0    0
[5,] 1.878e-10    0    0    0    0
```




How about at higher harvest levels?


```r
M[[4]][indices_fixed_x(5), indices_fixed_x(6)]
```



```
       [,1]    [,2] [,3]     [,4] [,5]
[1,] 0.2896 0.00000    0 0.000000    0
[2,] 0.2184 0.00000    0 0.000000    0
[3,] 0.0000 0.14756    0 0.000000    0
[4,] 0.0000 0.07719    0 0.000000    0
[5,] 0.0000 0.00000    0 0.007264    0
```



```r
M[[4]][indices_fixed_x(5), indices_fixed_x(10)]
```



```
          [,1] [,2] [,3] [,4] [,5]
[1,] 0.0255697    0    0    0    0
[2,] 0.0191809    0    0    0    0
[3,] 0.0128325    0    0    0    0
[4,] 0.0065240    0    0    0    0
[5,] 0.0002551    0    0    0    0
```






Static solution


```r
bevholt <- function(x,h, p) p[1] * (x-h) / (1 - p[2] * (x-h))
sdp <- determine_SDP_matrix(bevholt, c(1.5, 0.05), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
```





Confirm that the policy with high probability on model 1 matches the static solution for model 1:


```r
static$D
```



```
      [,1] [,2] [,3] [,4] [,5]
 [1,]    1    1    1    1    1
 [2,]    1    1    1    1    1
 [3,]    1    1    1    1    4
 [4,]    1    1    1    1    5
 [5,]    1    1    1    1    6
 [6,]    2    2    2    2    7
 [7,]    3    3    3    3    8
 [8,]    4    4    4    4    9
 [9,]    5    5    5    5   10
[10,]    6    6    6    6   10
```




Is the policy any different if most of our belief is on model 2?


```r
extract_policy(active$D, length(p_grid), length(x_grid), length(p_grid)) 
```



```
      [,1] [,2] [,3] [,4] [,5]
 [1,]    1    1    1    1    1
 [2,]    1    1    1    1    1
 [3,]    1    1    1    1    3
 [4,]    1    1    1    1    4
 [5,]    1    1    1    1    5
 [6,]    2    2    2    1    6
 [7,]    3    3    3    2    7
 [8,]    4    4    4    4    8
 [9,]    5    5    5    5    9
[10,]    6    6    6    6   10
```



```r
extract_policy(active$D, 1, length(x_grid), length(p_grid)) 
```



```
      [,1] [,2] [,3] [,4] [,5]
 [1,]    1    1    1    1    1
 [2,]    1    1    1    1    1
 [3,]    1    1    1    1    3
 [4,]    2    2    2    2    4
 [5,]    3    3    3    3    5
 [6,]    4    4    4    4    6
 [7,]    5    5    5    5    7
 [8,]    6    6    6    6    8
 [9,]    7    7    7    7    9
[10,]    8    8    8    8   10
```







```r
update_belief = function(x_t0, p_t0, x_t1){
  y1 = p_t0 * f1(x_t1, x_t0)
  y2 = (1-p_t0) * f2(x_t1, x_t0)
  P1 = y1 / (y1 + y2)
  if(is.na(P1) || x_t0 == 0)
    P1 = p_t0
  else{
    i = 1
    np = length(p_grid)
    while(p_grid[i] < P1 & i < np)
      i = i+1
    P1 = p_grid[i]  
 }
 P1
}


active_adaptive_simulate <- function(f, pars, x_grid, h_grid, p_grid, x0, 
                                     p0, D, z_g, update_belief){
  # initialize variables with initial conditions
  OptTime <- dim(D)[2]    # Stopping time
  x_h <- numeric(OptTime) # population dynamics with harvest
  h <- numeric(OptTime)   # optimal havest level
  x_h[1] <- x0            # initial values
  p <- numeric(OptTime)   # belief
  p[1] <- p0              # initial belief
  s <- x_h                # also track escapement
  x <- x_h                # What would happen with no havest
  nx <- length(x_grid)
  np <- length(p_grid)
  getpolicy <- function(p,x, time) D[x + (p-1)*nx, time]
  ## Simulate through time ##
  for(t in 1:(OptTime-1)){
    # Current state (is closest to which grid posititon) 
    St <- which.min(abs(x_grid - x_h[t])) 
    h[t] <- h_grid[getpolicy(p[t], St, t)] 
    # Implement harvest/(effort) based on quota with noise 
    # Noise in growth 
    z <- z_g() 
    # population grows
    x_h[t+1] <- z * f(x_h[t], h[t], pars) # with havest
    s[t]     <- x_h[t] - h[t] # anticipated escapement
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
    p[t+1]   <- update_belief(x[t+1], p[t], x[t]) 
  }
  # formats output 
  data.frame(time = 1:OptTime, fishstock = x_h, harvest = h,
             unharvested = x, escapement = s, belief = p) 
}
```






```r
z_g <- function() rlnorm(1,  0, sigma_g) 
sim <- active_adaptive_simulate(bevholt, c(1.5, 0.05), x_grid, h_grid, p_grid, 
                                x_grid[6], p_grid[4], active$D, z_g, update_belief)
require(ggplot2)
ggplot(sim) + geom_line(aes(time, fishstock)) + geom_line(aes(time, harvest), col="green") 
```

![plot of chunk unnamed-chunk-5](http://farm8.staticflickr.com/7209/6870442654_47e6fd987b_n.jpg) 

```r
ggplot(sim) + geom_line(aes(time, belief)) 
```

![plot of chunk unnamed-chunk-5](http://farm8.staticflickr.com/7217/7016551421_150373a9fa_n.jpg) 


