




 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

## Setup the system


```r
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```






```r
delta <- 0.05     # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 0       # bonus for satisfying the boundary condition
```








```r
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```








```r
f <- BevHolt                # Select the state equation
pars <- c(1.5, 0.05)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference 
xT <- 0                     # boundary conditions
x0 <- K
```






```r
profit <- profit_harvest(price = 10, c0 = 30, c1 = 10)
```






```r
x_grid <- seq(0.01, 1.2 * K, length = gridsize)  
h_grid <- seq(0.01, 0.8 * K, length = gridsize)  
```




## Declare the different penalty norms



```r
L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
asymmetric <- function(c2) function(h, h_prev)  c2 * max(h - h_prev, 0)
fixed <-  function(c2) function(h, h_prev) c2
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2
```




Calculate the transition matrix and intialize the the list of penalty functions and cost levels we will be looping over.


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
penaltyfns <- list(L2=L2, L1=L1, asy=asymmetric, fixed=fixed)
c2 <- seq(0, 30, length.out = 31)
```




This can take a while, so we use explicit parallelization, 


```r
require(snowfall)
sfInit(cpu=4, parallel=T)
```



```
R Version:  R version 2.14.2 (2012-02-29) 

```



```r
sfLibrary(pdgControl)
```



```
Library pdgControl loaded.
```



```r
sfExportAll()
```




## Loop over penalty functions and magnitudes



```r
policies <- 
sfSapply(penaltyfns, function(penalty){
  policies <- 
  sapply(c2, function(c2){
      policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                        profit, delta, reward, penalty = penalty(c2))
      i <- which(x_grid > K)[1]
      max(policycost$penalty_free_V[i,]) 
  })
})
```




Quadratic costs on fishing effort have to be done separately,


```r
quad <- 
  sapply(c2, function(c2){
  effort_penalty = function(x,h) .1*c2*h/x
  policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                        profit, delta, reward, penalty = fixed(0),     
                        effort_penalty)
      i <- which(x_grid > K)[1]
      max(policycost$penalty_free_V[i,]) # chooses the most sensible harvest in t=1
})
dat <- cbind(policies, quad)
```




Tidy up the data and plot the net present value (before the penalty has been paid) relative to that achieved when managed without a penalty.  


```r
npv0 <- dat[1,3] 
dat <- data.frame(c2=c2,dat)
dat <- melt(dat, id="c2")
ggplot(dat, aes(c2, (npv0-value)/npv0, col=variable)) + geom_point() + geom_line()
```

![plot of chunk unnamed-chunk-5](http://farm7.staticflickr.com/6233/6996636029_7e6a5a208d_o.png) 



