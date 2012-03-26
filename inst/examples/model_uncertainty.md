




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
sdp <- determine_SDP_matrix(Myers, c(1.1, 2, 10), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static_sim <- ForwardSimulate(Myers, c(1.1,2,10), x_grid, h_grid, 
                              K, static$D, z_g)
static$D
```



```
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
 [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [4,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [5,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [6,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [7,]    2    2    2    2    2    2    2    2    2     2     2     2     2
 [8,]    3    3    3    3    3    3    3    3    3     3     3     3     3
 [9,]    4    4    4    4    4    4    4    4    4     4     4     4     4
[10,]    6    6    6    6    6    6    6    6    6     6     6     6     6
[11,]    7    7    7    7    7    7    7    7    7     7     7     7     7
      [,14] [,15]
 [1,]     1     1
 [2,]     1     1
 [3,]     1     1
 [4,]     1     4
 [5,]     1     5
 [6,]     1     6
 [7,]     2     7
 [8,]     3     8
 [9,]     4    10
[10,]     5    11
[11,]     6    11
```



```r
static_sim
```



```
   time fishstock harvest unharvested escapement
1     1    10.000       6      10.000     4.0000
2     2     8.413       3      12.428     5.4126
3     3     8.342       3      10.509     5.3416
4     4     8.526       5      10.558     3.5256
5     5     6.405       1      10.607     5.4054
6     6     6.387       1       7.874     5.3875
7     7     7.415       2       8.585     5.4153
8     8     7.455       2       8.802     5.4545
9     9     7.126       2       8.433     5.1265
10   10     6.932       2       8.390     4.9317
11   11     7.483       2       9.247     5.4835
12   12     4.110       0       4.904     4.1103
13   13     8.388       3       9.431     5.3878
14   14    10.965      10      13.253     0.9655
15   15     1.107       0      12.281     0.0000
```




Active Adaptive Mangement solution


```r
bevholt <- function(x, h, p) max(p[1] * (x - h) / (1 - p[2] * (x - h)), 0)
myers  <- function(x, h, p) max(p[1] * (x - h) ^ 2 / (1 + (x - h) ^ 2 / p[2]), 0)
f1 <- setmodel(myers, c(1.1, 10))
f2 <- setmodel(bevholt, c(1.5, 0.05))
#f2 <- setmodel(bevholt, c(1.6, 0.05))

M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
```





What if we begin with great doubts about the allee model:


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

![plot of chunk activeplots](http://farm7.staticflickr.com/6092/6872772156_8d0ff58560_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7224/6872772468_f3798d7262_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm7.staticflickr.com/6111/7018879281_b2bc4ea6a2_o.png) 



Myers Static solution


```r
sdp <- determine_SDP_matrix(Myers, c(1.1, 2, 10), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
static_sim <- ForwardSimulate(Myers, c(1.1,2,10), x_grid, h_grid, 
                              K, static$D, z_g)
static$D
```



```
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
 [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [4,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [5,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [6,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [7,]    2    2    2    2    2    2    2    2    2     2     2     2     2
 [8,]    3    3    3    3    3    3    3    3    3     3     3     3     3
 [9,]    4    4    4    4    4    4    4    4    4     4     4     4     4
[10,]    6    6    6    6    6    6    6    6    6     6     6     6     6
[11,]    7    7    7    7    7    7    7    7    7     7     7     7     7
      [,14] [,15]
 [1,]     1     1
 [2,]     1     1
 [3,]     1     1
 [4,]     1     4
 [5,]     1     5
 [6,]     1     6
 [7,]     2     7
 [8,]     3     8
 [9,]     4    10
[10,]     5    11
[11,]     6    11
```



```r
static_sim
```



```
   time fishstock harvest unharvested escapement
1     1    10.000       6      10.000      4.000
2     2     5.680       1       8.390      4.680
3     3     6.334       1       8.079      5.334
4     4     7.011       2       8.216      5.011
5     5    10.374       6      12.634      4.374
6     6     7.957       3      11.402      4.957
7     7     6.607       2       8.633      4.607
8     8    11.262       6      14.608      5.262
9     9     6.597       2       8.577      4.597
10   10     6.275       1       8.138      5.275
11   11     6.943       2       8.200      4.943
12   12     6.403       1       7.856      5.403
13   13     8.490       3       9.809      5.490
14   14     7.567       9       9.129     -1.433
15   15     0.000       0       8.970      0.000
```




