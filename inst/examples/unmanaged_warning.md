




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
  a[t+1] = a[t] + .01
}
```






### Plot of timeseries 



```r
plot(x, type='l')
```

![plot of chunk p0](http://farm9.staticflickr.com/8307/7845763086_ff30bc53c4_o.png) 


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
m <- model$pars["m"]
```

```
Error: object 'model' not found
```





Set up a parallel environment





Evaluate the ROC curve



```r
clusterExport(cl, ls())
clusterExport(cl, list = c("A", "B"))
```

```
Error: object 'cl' not found
```

```r
reps <- parLapply(cl, 1:500, function(i) compare(A, B))
```

```
Error: object 'cl' not found
```

```r
lr <- lik_ratios(reps)
```

```
Error: object 'reps' not found
```

```r
roc <- roc_data(lr)
```

```
Error: object 'lr' not found
```




Plot results.



```r
require(ggplot2)
ggplot(lr) + geom_density(aes(value, fill = simulation), alpha = 0.6) + 
    geom_vline(aes(xintercept = observed))
```

```
Error: object 'lr' not found
```

```r
ggplot(roc) + geom_line(aes(False.positives, True.positives))
```

```
Error: object 'roc' not found
```





