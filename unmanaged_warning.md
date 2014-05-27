---
layout: page
---






# Detecting warning signals in an unmanaged system 

This is a precursor exercise to managed_warning.Rmd to calibrate a reasonable set of parameters over which to evaluate the detection.  






Simulate from the May population dynamics function





```r
n <- 5000
x <- vector(mode="double", length=n)
a <- x
x[1] <- 8 # positive equilibrium
z <- rlnorm(n, 0, .1)
r = .75; k = 10; a[1]=1; H=1; Q = 3
for(t in 1:n){
  x[t+1] = z[t] *  x[t] * exp(r * (1 - x[t] / k) - a[t] * x[t] ^ (Q - 1) / (x[t] ^ Q + H ^ Q)) 
  a[t+1] = a[t] + 1/5000
}
```






### Plot of timeseries 



```r
plot(x, type='l')
```

![plot of chunk p0](http://farm8.staticflickr.com/7108/7853908652_4ba7454256_o.png) 


Truncate the timeseries 



```r
y <- x[x > 1.5]
```






```r
plot(y, type='l')
```

![plot of chunk p1](http://farm9.staticflickr.com/8303/7853908906_b8a9307171_o.png) 


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
0.908 
```

```r
var_tau
```

```
    tau 
-0.5354 
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
[1] 119.8
```

```r
m
```

```
         m 
-2.658e-05 
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

![plot of chunk summarystat_roc](http://farm9.staticflickr.com/8441/7853909040_3a772c536f_o.png) 





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

![plot of chunk plotroc](http://farm8.staticflickr.com/7108/7853909166_45ecb1af08_o.png) 

```r
ggplot(roc) + geom_line(aes(False.positives, True.positives))
```

![plot of chunk plotroc](http://farm9.staticflickr.com/8448/7853909274_32a730a48a_o.png) 



