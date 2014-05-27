---
layout: page
---






System setup 


```r
require(pdgControl)
T <- 15
xT <- 0
reward <- 0
z_g <- function() rlnorm(1,  0, sigma_g) 
profit <- profit_harvest(price = 10, c0 = 1) 
delta <- 0.05
sigma_g = 0.05
m_pars <- c(1.5, 5)
K <-  .5 * (m_pars[1] * m_pars[2] + sqrt((m_pars[1] * m_pars[2]) ^ 2 - 4 * m_pars[2])) 
pars <- c(1.5, (1.5-1)/K)
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(.01,K,length=15) 
h_grid <- seq(0, K, length=11)
```




BevHolt Static solution


```r
sdp <- determine_SDP_matrix(BevHolt, pars, x_grid, h_grid, sigma_g)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static_sim <- ForwardSimulate(BevHolt, pars, x_grid, h_grid, 
                              K, static$D, z_g)
static$D
```



```
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
 [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [4,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [5,]    1    1    1    1    1    1    1    1    1     1     2     2     2
 [6,]    1    1    1    1    1    1    1    1    1     1     1     2     2
 [7,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [8,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [9,]    2    2    2    2    2    2    2    2    2     2     2     2     2
[10,]    3    3    3    3    3    3    3    3    3     3     3     3     3
[11,]    3    3    3    3    3    3    3    3    3     3     3     3     3
[12,]    5    5    5    5    5    5    5    5    5     5     5     5     5
[13,]    5    5    5    5    5    5    5    5    5     5     5     5     5
[14,]    5    5    5    5    5    5    5    5    5     5     5     5     5
[15,]    7    7    7    7    7    7    7    7    7     7     7     7     7
      [,14] [,15]
 [1,]     1     1
 [2,]     1     2
 [3,]     1     3
 [4,]     1     4
 [5,]     2     4
 [6,]     2     5
 [7,]     1     6
 [8,]     2     6
 [9,]     3     7
[10,]     3     8
[11,]     5     9
[12,]     5     9
[13,]     5    10
[14,]     5    11
[15,]     7    11
```



```r
static_sim
```



```
   time fishstock harvest unharvested escapement
1     1     6.760   4.056       6.760     2.7042
2     2     3.085   0.000       6.170     3.0848
3     3     3.802   0.676       6.413     3.1264
4     4     3.571   0.000       6.118     3.5713
5     5     3.861   0.676       5.757     3.1851
6     6     3.819   0.676       5.981     3.1425
7     7     3.915   0.676       6.366     3.2387
8     8     4.041   0.676       6.695     3.3654
9     9     4.359   1.352       7.244     3.0073
10   10     3.865   0.676       7.410     3.1890
11   11     3.890   0.676       7.216     3.2136
12   12     3.993   0.676       7.235     3.3168
13   13     3.632   1.352       6.427     2.2801
14   14     3.021   3.380       6.746    -0.3588
15   15     0.000   0.000       6.638     0.0000
```





Myers Static solution


```r
sdp <- determine_SDP_matrix(Myers, c(m_pars[1], 2, m_pars[2]), x_grid, h_grid, sigma_g)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static_sim <- ForwardSimulate(Myers,  c(m_pars[1], 2, m_pars[2]), x_grid, h_grid, 
                              K, static$D, z_g)
static$D
```



```
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
 [1,]    2    2    2    2    2    2    2    2    2     2     2     2     2
 [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [4,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [5,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [6,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [7,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [8,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [9,]    2    2    2    2    2    2    2    2    2     2     2     2     2
[10,]    3    3    3    3    3    3    3    3    3     3     3     3     3
[11,]    3    3    3    3    3    3    3    3    3     3     3     3     3
[12,]    4    4    4    4    4    4    4    4    4     4     4     4     4
[13,]    5    5    5    5    5    5    5    5    5     5     5     5     5
[14,]    5    5    5    5    5    5    5    5    5     5     5     5     5
[15,]    6    6    6    6    6    6    6    6    6     6     6     6     6
      [,14] [,15]
 [1,]     2     2
 [2,]     1     2
 [3,]     1     3
 [4,]     1     4
 [5,]     1     4
 [6,]     1     5
 [7,]     1     6
 [8,]     1     7
 [9,]     2     7
[10,]     3     8
[11,]     3     9
[12,]     4     9
[13,]     5    10
[14,]     5    11
[15,]     6    11
```



```r
static_sim
```



```
   time fishstock harvest unharvested escapement
1     1     6.760   3.380       6.760      3.380
2     2     5.023   1.352       6.509      3.671
3     3     5.653   2.704       6.932      2.949
4     4     4.933   1.352       7.038      3.581
5     5     5.652   2.704       7.135      2.948
6     6     4.606   1.352       6.607      3.254
7     7     4.733   1.352       6.252      3.381
8     8     5.091   2.028       6.488      3.063
9     9     4.815   1.352       6.598      3.463
10   10     5.642   2.704       7.170      2.937
11   11     4.721   1.352       6.796      3.369
12   12     4.992   1.352       6.489      3.640
13   13     5.471   2.028       6.736      3.443
14   14     5.563   6.084       7.125     -0.521
15   15     0.000   0.000       6.182      0.000
```






Active Adaptive Mangement solution


```r
bevholt <- function(x, h, p) max(p[1] * (x - h) / (1 - p[2] * (x - h)), 0)
myers  <- function(x, h, p) max(p[1] * (x - h) ^ 2 / (1 + (x - h) ^ 2 / p[2]), 0)
f1 <- setmodel(myers, m_pars)
f2 <- setmodel(bevholt, pars)

M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
```







```r
sims <- lapply(1:100, function(i){
  active_adaptive_simulate(Myers, c(m_pars[1], 2, m_pars[2]), x_grid, h_grid, p_grid, 
                                K, p_grid[5], active$D,
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

![plot of chunk activeplots](http://farm8.staticflickr.com/7106/7019133749_38309d90a2_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7231/7019134197_cb196dbdb1_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7070/7019134585_5c9f217b91_o.png) 



