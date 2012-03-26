




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
2     2    8.0121       3      11.836     5.0121
3     3    5.3367       0       6.964     5.3367
4     4    8.2498       3       9.241     5.2498
5     5    7.5536       3       9.215     4.5536
6     6    6.8976       2       9.147     4.8976
7     7    6.3917       1       8.090     5.3917
8     8    5.5196       1       6.435     4.5196
9     9    8.9194       5      10.702     3.9194
10   10    5.8773       1       8.924     4.8773
11   11    8.8304       5      11.143     3.8304
12   12    5.7452       1       8.941     4.7452
13   13    5.1690       0       6.635     5.1690
14   14    4.3829       4       4.908     0.3829
15   15    0.1168       0       5.715     0.0000
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

![plot of chunk activeplots](http://farm8.staticflickr.com/7243/6872765734_40fbbda7dd_o.png) 

```r

ggplot(dat) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm7.staticflickr.com/6114/6872766082_bcac37cbff_o.png) 

```r
ggplot(dat) + geom_line(aes(time, belief, group = reps), alpha = 0.2)
```

![plot of chunk activeplots](http://farm7.staticflickr.com/6045/6872766372_13d87b7db4_o.png) 



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
1     1    10.000       6      10.000     4.0000
2     2     6.721       2       9.928     4.7208
3     3     8.805       5      11.581     3.8050
4     4     5.249       0       8.259     5.2495
5     5     4.885       0       5.806     4.8848
6     6     6.386       1       6.989     5.3856
7     7     5.849       1       6.529     4.8494
8     8    10.189       6      11.762     4.1885
9     9     5.437       0       7.961     5.4370
10   10     7.503       3       8.673     4.5034
11   11     8.661       5      11.414     3.6609
12   12     6.533       2      10.595     4.5333
13   13     7.128       2       9.730     5.1284
14   14     8.117       9      10.133    -0.8832
15   15     0.000       0      15.811     0.0000
```




