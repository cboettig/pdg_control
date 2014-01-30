


# Reed Model

 * Author [Carl Boettiger](http://carlboettiger.info), <cboettig@gmail.com>
 * License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
 * Description:  Implements a numerical version of the SDP described in <a href="http://dx.doi.org/10.1016/0095-0696(79)90014-7">Reed (1979)</a> .





Chose the state equation / population dynamics function


```r
#f <- RickerAllee
#pars <- c(2, K, 5)
f <- BevHolt
pars <- c(2,4)
K <- (pars[1]-1)/pars[2]

#K <- 100
#pars <- c(1,K)
#f <- function(x,h,p){
#  sapply(x, function(x){
 #   S = max(x - h, 0)
#    p[1] * S * (1 - S/p[2]) + S
#  })
#}
```


We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  

<div> $$ \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x}, $$ </div> 


conditioned on h > x and x > 0,


```r
price <- 1
c0 <- 0
c1 <- 0
profit <- profit_harvest(price=price, c0 = c0, c1=c1) 
```


with price = 1, `c0` = 0 and `c1` = 0. 



```r
xmin <- 0
xmax <- K*1.5
grid_n <- 200
```


We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from 0 to 0.375 on a grid of 200 points, and over an identical discrete grid of possible harvest values.  



```r
x_grid <- seq(xmin, xmax, length = grid_n)  
h_grid <- x_grid  
```




```r
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- .1
```


We will determine the optimal solution over a 25 time step window with boundary condition for stock at 0 and discounting rate of 0.05.  The Reed model considers a stochastic growth model 

<div> $$ x_{t+1} = z_g f(x_t) $$ </div> 

for the random variable `z_g`, given by 


```r
#z_g <- function() 1 + rlnorm(1,0, sigma_g)  
pdfn <- function(P, s) dlnorm(P, 0, s)

#z_g <- function() 1+(2*runif(1, 0,  1)-1) * sigma_g
#pdfn <- function(P, s)  dunif(P, 1 - s, 1 + s)
```






```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g, pdfn)
```


### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward=0)
opt$S
```

```
[1] 0.1036
```




Stationary optimal policy:  


```r
s_opt <- value_iteration(SDP_Mat, x_grid, h_grid, OptTime=1000, xT, profit, delta)
```



## Compare to Clark (noise-free)


```r
SDP_Mat2 <- determine_SDP_matrix(f, pars, x_grid, h_grid, 0.001, pdfn)
```


Bellman's algorithm to compute the optimal solution for all possible trajectories.


```r
det <- find_dp_optim(SDP_Mat2, x_grid, h_grid, OptTime, xT, profit, delta, reward=0)
```




Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   


Plot the policy function (in terms of escapement, `x-h`, rather than harvest `h`) at equilibrium (first time-step):


```r
require(reshape2)
policies <- melt(data.frame(stock=x_grid, 
                            S = x_grid[opt$D[,1]], 
                            D = x_grid[det$D[,1]], 
                            Stationary = x_grid[s_opt$D]
                            ), id="stock")
q1 <- ggplot(policies, aes(stock, stock - value, color=variable)) + geom_point(alpha=.4) + xlab("stock size") + ylab("escapement") 
q1
```

![plot of chunk policyfn_plot](http://farm6.staticflickr.com/5512/12227370565_722ac24997_o.png) 


and the value function (at equilibrium):


```r
q2 <- qplot(x_grid, opt$V, xlab="stock size", ylab="value") + 
geom_vline(xintercept=opt$S)
q2
```

![plot of chunk valuefn_plot](http://farm3.staticflickr.com/2837/12227371395_1ff040d11b_o.png) 






