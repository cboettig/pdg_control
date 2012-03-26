




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
1     1     6.760   4.056       6.760    2.70416
2     2     3.327   0.000       6.655    3.32729
3     3     3.959   0.676       6.613    3.28332
4     4     4.039   0.676       6.789    3.36282
5     5     4.078   0.676       6.845    3.40225
6     6     3.981   0.676       6.655    3.30502
7     7     4.077   0.676       6.847    3.40093
8     8     4.085   0.676       6.833    3.40915
9     9     3.979   0.676       6.634    3.30301
10   10     3.737   0.676       6.265    3.06110
11   11     3.800   0.676       6.517    3.12357
12   12     3.988   0.676       6.912    3.31201
13   13     4.003   1.352       6.882    2.65063
14   14     3.290   3.380       6.770   -0.09044
15   15     0.000   0.000       6.706    0.00000
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
1     1     6.760   3.380       6.760     3.3802
2     2     5.065   1.352       6.563     3.7124
3     3     5.552   2.028       6.779     3.5235
4     4     5.401   2.028       6.833     3.3730
5     5     5.379   2.028       6.993     3.3505
6     6     5.204   2.028       6.824     3.1757
7     7     4.955   1.352       6.693     3.6028
8     8     5.502   2.028       6.856     3.4737
9     9     5.035   1.352       6.436     3.6825
10   10     5.919   2.704       7.229     3.2146
11   11     5.171   2.028       7.004     3.1433
12   12     4.903   1.352       6.701     3.5510
13   13     5.218   2.028       6.557     3.1902
14   14     4.975   5.408       6.646    -0.4333
15   15     0.000   0.000       7.084     0.0000
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

![plot of chunk activeplots](http://farm8.staticflickr.com/7062/7019126121_b597fa14a8_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7046/7019126361_6f14662ffe_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7227/7019126511_9fb725d60a_o.png) 



