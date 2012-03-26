




System setup 


```r
require(pdgControl)
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(1,10,length=10) 
sigma_g = 0.2
h_grid <- seq(0, 10, length=11)
T <- 15
xT <- 0
z_g <- function() rlnorm(1,  0, sigma_g) 
profit <- profit_harvest(price = 10, c0 = 30) 
delta <- 0.05
reward <- 0
pars <- c(1.5, 0.05)
K <- (pars[1]-1)/pars[2]
```




Static solution


```r
sdp <- determine_SDP_matrix(BevHolt, pars, x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static_sim <- ForwardSimulate(BevHolt, pars, x_grid, h_grid, 
                              K, static$D, z_g)
```




Active Adaptive Mangement solution


```r
bevholt <- function(x, h, p) max(p[1] * (x - h) / (1 - p[2] * (x - h)), 0)
myers  <- function(x, h, p) max(p[1] * (x - h) ^ 2 / (1 - (x - h) ^ 2 / p[2]), 0)
#f1 <- setmodel(myers, c(1.5, 10))
f1 <- setmodel(bevholt, c(1.5, 0.05))
f2 <- setmodel(bevholt, c(4, 0.05))

M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
```







```r
sim <- active_adaptive_simulate(BevHolt, pars, x_grid, h_grid, p_grid, 
                                K, p_grid[1], active$D,
                                z_g, update_belief(f1,f2))
require(ggplot2)
ggplot(sim) + geom_line(aes(time, fishstock)) + geom_line(aes(time, harvest), col="green") 
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7067/7018496791_a1c4687b3d_o.png) 

```r
ggplot(sim) + geom_line(aes(time, belief)) 
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7037/6872389354_47c5d1eeb8_o.png) 




