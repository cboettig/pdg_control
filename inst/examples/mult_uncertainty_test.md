Test multiple uncertainty function

  

```
## Loading required package: pdgControl
```

```
## Loading required package: reshape2
```

```
## Loading required package: ggplot2
```

```
## Loading required package: data.table
```






## Current version of the multiple uncertainty function


```r
#' SDP under multiple uncertainty
#' 
#' Computes the SDP solution under the case of growth noise,
#' implementation errors in harvesting, and meaurement errors in the stock
#' assessment.
#' 
#' @param f the growth function of the escapement population (x-h) should
#' be a function of f(t, y, p), with parameters p @param p the parameters
#' of the growth function @param x_grid the discrete values allowed for
#' the population size, x @param h_grid the discrete values of harvest
#' levels to optimize over @param sigma is the shape parameters for noise
#' distribution (sigma_g, sigma_m, sigma_i) (default is no noise) @param
#' pdfn is the probability density function (same functional form is used
#' for growth, measure, implement).  (Default is uniform) @param profit is
#' the profit function (defaults to the realized harvest) @return The D
#' matrix giving the optimal harvest for each possible state in each
#' timestep.  Hence the optimal policy at time t is given by the policy
#' function D[,t].  Values of D[,t] correspond to the index of h_grid.
#' Indices of of D[,t] correspond to states in y_grid.  @export
SDP_multiple_uncertainty <- function(f, p, x_grid, h_grid, Tmax = 25, 
    sigmas = c(sigma_g = 0, sigma_m = 0, sigma_i = 0), pdfn = pdfn, profit = function(x, 
        h) pmin(x, h)) {
    
    ## Initialize stuff
    sigma_g <- sigmas["sigma_g"]
    sigma_m <- sigmas["sigma_m"]
    sigma_i <- sigmas["sigma_i"]
    n_x <- length(x_grid)
    n_h <- length(h_grid)
    D <- matrix(NA, nrow = n_x, ncol = Tmax)
    
    ## Compute the transition matrices
    P <- outer(x_grid, h_grid, profit)
    M <- rownorm(outer(x_grid, x_grid, pdfn, sigma_m))
    I <- rownorm(outer(h_grid, h_grid, pdfn, sigma_i))
    
    # integrate over uncertainty in I and X
    F <- lapply(1:n_h, function(q) {
        t(sapply(1:n_x, function(y) {
            out <- numeric(n_x)
            mu <- sum(sapply(1:n_x, function(x) f(x_grid[x], h_grid, p) %*% 
                I[q, ] * M[x, y]))  # Implementation & Measurement error
            if (mu == 0) 
                out[1] <- 1 else out <- dlnorm(x_grid/mu, 0, sigma_g)
            out/sum(out)
        }))
    })
    
    Ep <- M %*% P %*% I
    V <- Ep
    for (t in 1:Tmax) {
        D[, (Tmax - t + 1)] <- apply(V, 1, which.max)
        v_t <- apply(V, 1, max)  # vector of current values
        V <- sapply(1:n_h, function(j) {
            # updated value matrix
            Ep[, j] + (1 - delta) * M %*% F[[j]] %*% v_t
        })
    }
    list(D = D, M = M, I = I, P = P, Ep = Ep, V = V, F = F)
}

# row-normalize the probability distribution
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))
# binning routine (vectorized)
snap_to_grid <- function(x, grid) sapply(x, function(x) grid[which.min(abs(grid - 
    x))])
```


 


```r
f <- function(x, h, p) {
    A <- p[1]
    B <- p[2]
    s <- pmax(x - h, 0)
    A * s/(1 + B * s)
}
profit = function(x, h) pmin(x, h)


pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
xmin <- 0
xmax <- 1.5 * K
n_x <- 120
n_h <- n_x
x_grid <- seq(xmin, xmax, length = n_x)
h_grid <- seq(xmin, xmax, length = n_h)
delta <- 0.05
xT <- 0
OptTime <- 25
sigma_g <- 0.5
```



  
Old school way:


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward = 0)
```



New way: 


```r
# Uniform pdf where width scales with mean (probably better way to do
# that) And handles degenerate/delta fn cases such as no width and zero
# mean
FUN <- function(P, mu, s) {
    if (mu == 0) {
        as.integer(P == 0)
    } else if (s > 0) {
        if (mu > 0) 
            dunif(P, mu * (1 - s), mu * (1 + s))
    } else {
        # delta spike
        P <- snap_to_grid(P, x_grid)
        mu <- snap_to_grid(mu, x_grid)
        as.integer(P == mu)
    }
}
pdfn <- Vectorize(FUN)
```



```r
new <- SDP_multiple_uncertainty(f, pars, x_grid, h_grid, OptTime, 
    sigmas = c(sigma_g = sigma_g, sigma_m = 0, sigma_i = 0), pdfn = pdfn)
m <- SDP_multiple_uncertainty(f, pars, x_grid, h_grid, OptTime, sigmas = c(sigma_g = sigma_g, 
    sigma_m = 0.1, sigma_i = 0), pdfn = pdfn)
i <- SDP_multiple_uncertainty(f, pars, x_grid, h_grid, OptTime, sigmas = c(sigma_g = sigma_g, 
    sigma_m = 0, sigma_i = 0.1), pdfn = pdfn)
```


Plot the policy function (in terms of escapement, `x-h`, rather than harvest `h`) at equilibrium (first time-step):


```r
require(reshape2)
policies <- melt(data.frame(stock = x_grid, old = x_grid[opt$D[, 
    1]], new = x_grid[new$D[, 1]], meas = x_grid[m$D[, 1]], imp = x_grid[i$D[, 
    1]]), id = "stock")
```



```r
q1 <- ggplot(policies, aes(stock, stock - value, color = variable)) + 
    geom_point() + xlab("stock size") + ylab("escapement")
q1
```

![plot of chunk policyfunctions](http://carlboettiger.info/assets/figures/2012-11-08-1c68e48d9a-policyfunctions.png) 



```r
rbind(opt$D[, 1], new$D[, 1], m$D[, 1], i$D[, 1])
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
## [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1
## [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1
## [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1
## [4,]    1    1    1    1    1    1    1    1    1     1     1     1     1
##      [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
## [1,]     1     1     1     1     1     1     1     1     1     1     1
## [2,]     1     1     1     1     1     1     1     1     1     1     1
## [3,]     1     1     1     1     1     1     1     1     1     1     1
## [4,]     1     1     1     1     1     1     1     1     1     1     1
##      [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
## [1,]     1     1     1     1     1     1     1     1     1     1     1
## [2,]     1     1     1     1     1     1     1     1     1     2     3
## [3,]     1     1     1     1     1     1     1     1     1     2     3
## [4,]     1     1     1     1     1     1     1     1     1    11    11
##      [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46]
## [1,]     2     3     4     5     6     7     8     9    10    11    12
## [2,]     4     5     6     7     8     9    10    11    12    13    14
## [3,]     4     6     7     8     9    10    11    12    12    13    14
## [4,]    11    11    11    11    11    11    11    11    11    11    11
##      [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55] [,56] [,57]
## [1,]    13    14    15    16    17    18    19    20    21    22    23
## [2,]    15    16    17    18    19    20    21    22    23    24    25
## [3,]    16    17    18    19    20    21    22    23    24    25    26
## [4,]    11    11    11    21    21    21    21    21    21    21    21
##      [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66] [,67] [,68]
## [1,]    24    25    26    27    28    29    30    31    32    33    34
## [2,]    26    27    28    29    30    31    32    33    34    35    36
## [3,]    27    28    29    30    31    32    33    34    35    35    36
## [4,]    21    21    21    31    31    31    31    31    31    31    31
##      [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77] [,78] [,79]
## [1,]    35    36    37    38    39    40    41    42    43    44    45
## [2,]    37    38    39    40    41    42    43    44    45    46    47
## [3,]    37    39    40    41    42    43    44    45    46    46    47
## [4,]    31    31    41    41    41    41    41    41    41    41    41
##      [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90]
## [1,]    46    47    48    49    50    51    52    53    54    55    56
## [2,]    48    49    50    51    52    53    54    55    56    57    58
## [3,]    48    49    51    52    53    54    55    56    57    57    58
## [4,]    41    41    41    41    41    41    41    41    52    61    61
##      [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100] [,101]
## [1,]    57    58    59    60    61    62    63    64    65     66     67
## [2,]    59    60    61    62    63    64    65    66    67     68     69
## [3,]    59    61    62    63    64    65    66    67    68     69     70
## [4,]    61    61    61    61    61    61    61    61    61     61     61
##      [,102] [,103] [,104] [,105] [,106] [,107] [,108] [,109] [,110] [,111]
## [1,]     68     69     70     71     72     73     74     75     76     77
## [2,]     70     71     72     73     74     75     76     77     78     79
## [3,]     71     72     73     73     74     75     76     77     77     77
## [4,]     61     61     61     61     72     72     81     81     81     81
##      [,112] [,113] [,114] [,115] [,116] [,117] [,118] [,119] [,120]
## [1,]     78     79     80     81     82     83     84     85     86
## [2,]     80     81     82     83     84     85     86     87     88
## [3,]     78     78     79     80     80     81     81     82     82
## [4,]     81     81     81     81     81     81     81    110    110
```

