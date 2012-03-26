




System setup 


```r
require(pdgControl)
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(.01,2.5,length=15) 
sigma_g = 0.2
h_grid <- seq(0, 2, length=11)
T <- 15
xT <- 0
z_g <- function() rlnorm(1,  0, sigma_g) 
profit <- profit_harvest(price = 10, c0 = 1) 
delta <- 0.05
reward <- 0
pars <- c(1.5, 0.25)
K <- (pars[1]-1)/pars[2]
```




BevHolt Static solution


```r
sdp <- determine_SDP_matrix(BevHolt, pars, x_grid, h_grid, .2)
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
 [5,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [6,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [7,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [8,]    2    2    2    2    2    2    2    2    2     2     2     2     2
 [9,]    3    3    3    3    3    3    3    3    3     3     3     3     3
[10,]    4    4    4    4    4    4    4    4    4     4     4     4     4
[11,]    5    5    5    5    5    5    5    5    5     5     5     5     5
[12,]    6    6    6    6    6    6    6    6    6     6     6     6     6
[13,]    7    7    7    7    7    7    7    7    7     7     7     7     7
[14,]    8    8    8    8    8    8    8    8    8     8     8     8     8
[15,]    9    9    9    9    9    9    9    9    9     9     9     9     9
      [,14] [,15]
 [1,]     1     1
 [2,]     1     2
 [3,]     1     3
 [4,]     1     4
 [5,]     1     5
 [6,]     1     6
 [7,]     1     7
 [8,]     2     8
 [9,]     3     9
[10,]     4     9
[11,]     5    10
[12,]     6    11
[13,]     7    11
[14,]     7    11
[15,]     8    11
```



```r
static_sim
```



```
   time fishstock harvest unharvested escapement
1     1     2.000     1.0       2.000     1.0000
2     2     0.906     0.0       1.510     0.9060
3     3     1.216     0.2       1.805     1.0162
4     4     1.637     0.6       2.512     1.0370
5     5     1.116     0.0       2.092     1.1162
6     6     1.402     0.4       2.206     1.0017
7     7     1.366     0.4       2.424     0.9657
8     8     1.067     0.0       2.071     1.0674
9     9     1.576     0.6       2.552     0.9758
10   10     1.100     0.0       2.185     1.0999
11   11     1.007     0.0       1.650     1.0073
12   12     1.008     0.0       1.463     1.0079
13   13     1.577     0.6       2.098     0.9765
14   14     1.020     1.2       1.788    -0.1802
15   15     0.000     0.0       2.315     0.0000
```







Active Adaptive Mangement solution


```r
bevholt <- function(x, h, p) max(p[1] * (x - h) / (1 - p[2] * (x - h)), 0)
myers  <- function(x, h, p) max(p[1] * (x - h) ^ 2 / (1 + (x - h) ^ 2 / p[2]), 0)
f1 <- setmodel(myers, c(1.1, 2))
f2 <- setmodel(bevholt, pars)

M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 
```







```r
sims <- lapply(1:100, function(i){
  active_adaptive_simulate(Myers, c(1.1, 2, 2), x_grid, h_grid, p_grid, 
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

![plot of chunk activeplots](http://farm8.staticflickr.com/7045/7019024915_481f321030_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm7.staticflickr.com/6117/6872919802_3a75aa1e07_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7081/7019025483_4a2689313f_o.png) 



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
 [2,]    2    2    2    2    2    2    2    2    2     2     2     2     2
 [3,]    3    3    3    3    3    3    3    3    3     3     3     3     3
 [4,]    4    4    4    4    4    4    4    4    4     4     4     4     4
 [5,]    5    5    5    5    5    5    5    5    5     5     5     5     5
 [6,]    1    1    1    1    1    1    1    1    1     1     1     6     6
 [7,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [8,]    1    1    1    1    1    1    1    1    1     1     1     1     1
 [9,]    1    1    1    1    1    1    1    1    1     1     1     1     1
[10,]    1    1    1    1    1    1    1    1    1     1     1     1     1
[11,]    2    2    2    2    2    2    2    2    2     2     2     2     2
[12,]    3    3    3    3    3    3    3    3    3     3     3     3     3
[13,]    3    3    3    3    3    4    4    4    4     4     4     4     4
[14,]    4    4    4    4    4    4    4    4    4     4     4     4     4
[15,]    5    5    5    5    5    5    5    5    5     5     5     5     5
      [,14] [,15]
 [1,]     1     1
 [2,]     2     2
 [3,]     3     3
 [4,]     4     4
 [5,]     5     5
 [6,]     6     6
 [7,]     1     7
 [8,]     1     8
 [9,]     1     9
[10,]     2     9
[11,]     2    10
[12,]     3    11
[13,]     4    11
[14,]     5    11
[15,]     6    11
```



```r
static_sim
```



```
   time fishstock harvest unharvested escapement
1     1     2.000     0.4       2.000      1.600
2     2     1.504     0.0       2.108      1.504
3     3     1.189     0.0       1.983      1.189
4     4     1.227     0.0       2.797      1.227
5     5     1.259     0.0       4.223      1.259
6     6     1.863     0.2       8.721      1.663
7     7     2.404     0.6       9.811      1.804
8     8     2.760     0.8      10.181      1.960
9     9     3.957     0.8      13.004      3.157
10   10     5.992     0.8      11.335      5.192
11   11     8.698     0.8      11.064      7.898
12   12    10.108     0.8      10.842      9.308
13   13    16.475     1.0      16.936     15.475
14   14     7.830     2.0       7.882      5.830
15   15     9.234     0.0      10.294      0.000
```




