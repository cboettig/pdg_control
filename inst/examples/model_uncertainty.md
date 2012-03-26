




System setup 


```r
require(pdgControl)
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(.01,10,length=11) 
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
f1 <- setmodel(myers, c(1.1, 10))
f2 <- setmodel(bevholt, c(1.5, 0.05))
#f2 <- setmodel(bevholt, c(1.6, 0.05))

M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
```




Begin with minimal belief in the true model. 


```r
sims <- lapply(1:100, function(i){
  active_adaptive_simulate(Myers, c(1.1, 2, 10), x_grid, h_grid, p_grid, 
                                K, p_grid[1], active$D,
                                z_g, update_belief(f1,f2))
})
require(reshape2)
dat <- melt(sims, id=names(sims[[1]])) 
names(dat)[7] <- "reps"
require(ggplot2)
ggplot(subset(dat,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_line(aes(time, harvest), col="darkgreen") +  
  geom_line(aes(time, belief), col="darkred")
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7081/6872719126_d1df0ffd66_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm7.staticflickr.com/6031/7018825817_43381759c6_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7040/7018826047_263177a66f_o.png) 



