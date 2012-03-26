




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
1     1   10.0000       6      10.000     4.0000
2     2    7.0602       2      10.430     5.0602
3     3    6.3427       1       8.077     5.3427
4     4    6.9231       2       8.106     4.9231
5     5    7.0455       2       8.638     5.0455
6     6    9.7909       6      12.025     3.7909
7     7    5.5718       1       8.838     4.5718
8     8    5.3161       0       6.967     5.3161
9     9    7.6339       3       8.570     4.6339
10   10    7.3682       2       9.505     5.3682
11   11   10.0920       6      12.239     4.0920
12   12    5.9156       1       8.857     4.9156
13   13    5.1378       0       6.443     5.1378
14   14    7.4827       7       8.314     0.4827
15   15    0.2806       0      10.764     0.0000
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






```r
sims <- lapply(1:100, function(i){
  active_adaptive_simulate(Myers, c(1.1, 2, 10), x_grid, h_grid, p_grid, 
                                K, p_grid[4], active$D,
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

![plot of chunk activeplots](http://farm7.staticflickr.com/6036/7018850251_9b49ef7a8b_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm8.staticflickr.com/7102/7018850545_5562557904_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm7.staticflickr.com/6113/6872744460_b1976cd8f6_o.png) 



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
2     2     8.453       3      12.488      5.453
3     3    10.666       6      13.393      4.666
4     4     7.773       3      10.745      4.773
5     5     6.808       2       9.015      4.808
6     6     7.193       2       9.176      5.193
7     7     7.728       3       9.468      4.728
8     8     9.327       5      12.146      4.327
9     9     6.642       2       9.542      4.642
10   10     8.163       3      10.769      5.163
11   11     7.120       2       9.014      5.120
12   12     7.204       2       8.862      5.204
13   13     7.151       2       8.685      5.151
14   14     8.759      10      10.648     -1.241
15   15     0.000       0      15.285      0.000
```




