




System setup 


```r
require(pdgControl)
p_grid = seq(0.01,.99, length=10) 
x_grid = seq(1,10,length=15) 
sigma_g = 0.2
h_grid <- x_grid-1 
T <- 15
xT <- 0
z_g <- function() rlnorm(1,  0, sigma_g) 
profit <- profit_harvest(price=10, c0=30) 
delta <- 0.05
reward <- 0
bevholt <- function(x,h, p) p[1] * (x - h) / (1 - p[2] * (x - h))
myers  <- function(x,h,p) p[1] * (x - h) ^ 2 / (1 - (x - h) ^ 2 / p[2])
```




Static solution


```r
sdp <- determine_SDP_matrix(myers, c(1.5, 10), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
```




Active Adaptive Mangement solution


```r
f1 <- setmodel(myers, c(1.5,10))
f2 <- setmodel(bevholt, c(1.5,0.05))

M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
sim <- active_adaptive_simulate(myers, c(1.5, 10), x_grid, h_grid, p_grid, 
                                x_grid[6], p_grid[4], active$D, z_g, update_belief(f1,f2))
require(ggplot2)
ggplot(sim) + geom_line(aes(time, fishstock)) + geom_line(aes(time, harvest), col="green") 
```

![plot of chunk active](http://farm8.staticflickr.com/7259/6871560134_3b82478c08_o.png) 

```r
ggplot(sim) + geom_line(aes(time, belief)) 
```



```
Error: object 'belief' not found
```






