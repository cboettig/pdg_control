---
layout: page
---

---
layout: page
---





```
Error: .onLoad failed in loadNamespace() for 'XML', details: call:
dyn.load(file, DLLpath = DLLpath, ...)  error: unable to load shared
object
'/home/cboettig/R/x86_64-redhat-linux-gnu-library/2.15/XML/libs/XML.so':
libxmlsec1.so.1: cannot open shared object file: No such file or directory
```

```
Error: .onLoad failed in loadNamespace() for 'XML', details: call:
dyn.load(file, DLLpath = DLLpath, ...)  error: unable to load shared
object
'/home/cboettig/R/x86_64-redhat-linux-gnu-library/2.15/XML/libs/XML.so':
libxmlsec1.so.1: cannot open shared object file: No such file or directory
```




# Critical transition 

 * Author [Carl Boettiger](http://carlboettiger.info), <cboettig@gmail.com>
 * License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
 * Description:   





Chose the state equation / population dynamics function



```r
f <- May
```




With parameters 



```r
pars <- c(r = .75, k = 10, a=1, H=1, Q = 3)
K <- 8 # approx
```




Note that this bifurcates when a increases to 2.  


We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  

<div> $$ \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x}, $$ </div> 

conditioned on h > x and x > 0,



```r
price <- 1
c0 <- 0.0
c1 <- 0
profit <- profit_harvest(price=price, c0 = c0, c1=c1) 
```




with price = `1`, `c0` = `0` and `c1` = `0`. 




```r
xmin <- 0
xmax <- 1.5 * K
grid_n <- 100
```




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `12` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




```r
x_grid <- seq(xmin, xmax, length = grid_n)  
h_grid <- x_grid  
```







```r
delta <- 0.005
xT <- 0
OptTime <- 5000
sigma_g <- .1
```




We will determine the optimal solution over a `5000` time step window with boundary condition for stock at `0` and discounting rate of `0.005`.  The Reed model considers a stochastic growth model 

<div> $$ x_{t+1} = z_g f(x_t) $$ </div> 

for the random variable `z_g`, given by 



```r
z_g <- function() rlnorm(1, 0, sigma_g) 
```




No other sources of noise enter into the dynamics.  



```r
z_m <- function() 1
z_i <- function() 1
```








```r
pdfn <- function(P, s){
  dunif(P, 1 - s, 1 + s)
}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g, pdfn)
```





### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories.



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=OptTime, xT=xT, 
                     profit, delta=delta, reward=0)
```




Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   


Plot the policy function (in terms of escapement, `x-h`, rather than harvest `h`) at equilibrium (first time-step):



```r
q1 <- qplot(x_grid, x_grid - x_grid[opt$D[,1]], xlab="stock size", ylab="escapement") + 
geom_point(aes(x,y), data=data.frame(x=opt$S, y=opt$S), col="red")
q1
```

![plot of chunk policyfn_plot](http://farm9.staticflickr.com/8316/7889663888_6c2cc2de08_o.png) 


and the value function (at equilibrium):



```r
q2 <- qplot(x_grid, opt$V, xlab="stock size", ylab="value") + 
geom_vline(xintercept=opt$S)
q2
```

![plot of chunk valuefn_plot](http://farm9.staticflickr.com/8309/7889664294_939223f8ac_o.png) 


### Simulate 




```r
Dt <- 1

ForwardSimulate <- 
  function(f, pars, x_grid, h_grid, x0, D, z_g,
         z_m=function(x) 1, z_i = function(x) 1, 
         profit=NULL){

  OptTime <- dim(D)[2]    # Stopping time
  x_h <- numeric(OptTime) # population dynamics with harvest
  h <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0  # initial values
  s <- x_h # also track escapement
  x <- x_h # What would happen with no havest
  p <- numeric(OptTime)

  ## Simulate through time ##
  for(t in 1:(OptTime-1)){

    # Move towards collapse
    pars[3] <- pars[3] + Dt/OptTime



    # Assess stock, with potential measurement error
    m_t <- x_h[t] * z_m()
    # Current state (is closest to which grid posititon) 
    St <- which.min(abs(x_grid - m_t)) 
    q_t <- h_grid[D[St, (t + 1) ]] 
    # Implement harvest/(effort) based on quota with noise 
    h[t] <- q_t * z_i()
    # Noise in growth 
    z <- z_g() 
    # population grows
    x_h[t+1] <- z * f(x_h[t], h[t], pars) # with havest
    s[t+1]   <- x_h[t] - q_t # anticipated escapement
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
    p[t] <- profit(x_h[t], h[t])
  }
  data.frame(time = 1:OptTime, fishstock = x_h, harvest = h,
             unharvested = x, escapement = s, profit = p) 
}

dat <- ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i)
x <- dat$fishstock
```






### Plot of timeseries 



```r
plot(x, type='l')
```

![plot of chunk p0](http://farm9.staticflickr.com/8451/7889664498_8da4534b28_o.png) 


Truncate the timeseries 



```r
y <- x[x > 1.5]
```






```r
plot(y, type='l')
```

![plot of chunk p1](http://farm9.staticflickr.com/8177/7889664738_ae6dc140d5_o.png) 


### Calculate warning signals on the truncated series. 



```r
library(earlywarning)
dat <- data.frame(time=1:length(y), value=y)
```






```r
acor_tau <- warningtrend(dat, window_autocorr)
var_tau <- warningtrend(dat, window_var)

acor_tau
```

```
   tau 
0.9447 
```

```r
var_tau
```

```
   tau 
0.6166 
```




## Model based statistics


Fit the models



```r
A <- stability_model(dat, "OU")
B <- stability_model(dat, "LSN")
observed <- -2 * (logLik(A) - logLik(B))
m <- B$pars["m"]


observed
```

```
[1] 344
```

```r
m
```

```
         m 
-0.0001701 
```






Compute summary stat versions



```r
summarystat_roc <- function(A,B, summarystat_functions, reps=200){
  require(plyr)
  require(reshape2)
  Asim <- simulate(A, reps)
  Bsim <- simulate(B, reps)
  Asim <- melt(Asim, id = "time")
  Bsim <- melt(Bsim, id = "time")
  names(Asim)[2] <- "rep"
  names(Bsim)[2] <- "rep"
	dat <- lapply(summarystat_functions, function(f){
	  wsA <- ddply(Asim, "rep", warningtrend, f)
  	wsB <- ddply(Bsim, "rep", warningtrend, f)
	  data.frame(null = wsA$tau, test = wsB$tau)
	})
	tidy <- melt(dat)
}

dat <- summarystat_roc(A,B, list(var=window_var, acor=window_autocorr))
ggplot(dat) + geom_density(aes(value, fill=variable)) + facet_wrap(~L1)
```

![plot of chunk summarystat_roc](http://farm9.staticflickr.com/8035/7889664938_34fa3fff95_o.png) 





Set up a parallel environment



```r
require(snowfall)
sfInit(par=T, cpu=12)
```

```
R Version:  R version 2.15.0 (2012-03-30) 

```

```r
sfLibrary(earlywarning)
```

```
Library earlywarning loaded.
```

```r
sfExportAll() 
```




Evaluate the ROC curve



```r
reps <- sfLapply(1:500, function(i) compare(A, B))
lr <- lik_ratios(reps)
roc <- roc_data(lr)
```




Plot results.



```r
require(ggplot2)
ggplot(lr) + geom_density(aes(value, fill = simulation), alpha = 0.6) + 
    geom_vline(aes(xintercept = observed))
```

![plot of chunk plotroc](http://farm9.staticflickr.com/8445/7889665150_d2f03f6f86_o.png) 

```r
ggplot(roc) + geom_line(aes(False.positives, True.positives))
```

![plot of chunk plotroc](http://farm9.staticflickr.com/8442/7889665336_993d71c5d1_o.png) 




