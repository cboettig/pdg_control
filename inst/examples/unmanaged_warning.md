




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

![plot of chunk p0](http://farm9.staticflickr.com/8301/7847417972_d9a11d417d_o.png) 


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






Set up a parallel environment



```r
library(snow)
library(methods)
cl <- makeCluster(20, type = "MPI")
```

```
	20 slaves are spawned successfully. 0 failed.
```

```r
clusterEvalQ(cl, library(earlywarning))
```

```
[[1]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[2]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[3]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[4]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[5]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[6]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[7]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[8]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[9]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[10]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[11]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[12]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[13]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[14]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[15]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[16]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[17]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[18]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[19]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

[[20]]
 [1] "earlywarning" "snow"         "Rmpi"         "methods"     
 [5] "stats"        "graphics"     "grDevices"    "utils"       
 [9] "datasets"     "base"        

```






