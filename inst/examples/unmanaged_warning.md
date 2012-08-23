




# Detecting warning signals in an unmanaged system 

This is a precursor exercise to managed_warning.Rmd to calibrate a reasonable set of parameters over which to evaluate the detection.  






Simulate from the May population dynamics function





```r
n <- 5000
x <- vector(mode="double", length=n)
a <- x
x[1] <- 8 # positive equilibrium
z <- rlnorm(n, 0, .1)
r = .75; k = 10; a[1]=1.55; H=1; Q = 3
for(t in 1:n){
  x[t+1] = z[t] *  x[t] * exp(r * (1 - x[t] / k) - a[t] * x[t] ^ (Q - 1) / (x[t] ^ Q + H ^ Q)) 
  a[t+1] = a[t] + .001
}
```






### Plot of timeseries 



```r
plot(x, type='l')
```

![plot of chunk p0](http://farm9.staticflickr.com/8426/7847325690_d646752cb4_o.png) 


Truncate the timeseries 



```r
y <- x[x > 1.5]
```




### Calculate warning signals on the truncated series. 



```r
library(earlywarning)
dat <- data.frame(time=1:length(y), value=y)
```






```r
acor_tau <- warningtrend(dat, window_autocorr)
var_tau <- warningtrend(dat, window_var)
```




## Model based statistics


Fit the models



```r
A <- stability_model(dat, "OU")
B <- stability_model(dat, "LSN")
observed <- -2 * (logLik(A) - logLik(B))
m <- B$pars["m"]
```






