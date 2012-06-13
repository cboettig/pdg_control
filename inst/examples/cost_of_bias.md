






# Calculating the cost of bias  
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0


 * knitr-formatted [source code](https://github.com/cboettig/pdg_control/blob/master/inst/examples/cost_of_bias.Rmd)
 * [Cached data](http://two.ucdavis.edu/cboettig/data/cost_of_bias/)

Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).  Then compute the optimal solution under different forms of uncertainty and compare the results.  





## Model setup 

We will assume a Beverton-Holt state equation / population dynamics function, <span> \( f(x,h) = \frac{A x}{1 + B x} \)</span>



with parameters A = `1.5` and B = `0.05`.  The positive stationary root of the model is given by <span>\( \frac{A-1}{B} \)</span>, `10`.   




We also assume a profit function of the form <span>\( \Pi = p h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \)</span>, conditioned on <span>\( h > x \)</span> and <span>\(x > 0 \)</span>, with price p = `1`, c0 = `0.01`, and c1 = `0`.  





and solve the problem on a discrete grid of `100` for stock size and range `0`, `20`.  We use the same set of gridpoints for the possible harvest levels. 


## Scenarios 

We calculate the stochastic transition matrix for the probability of going from any state \\(x_t \\) to any other state \\(x_{t+1}\\) the following year, for each possible choice of harvest \\( h_t \\).  This provides a look-up table for the dynamic programming calculations.

### No Uncertainty 

The first scenario considers the completely deterministic case.  










We simulate 100 replicates of this system.  We will used a fixed seed so that we can compare these replicates to simulations under different conditions.  (Of course the seed is irrelevant at this stage since this is actually deterministic).  



### Growth uncertainty 




The next scenario introduces growth uncertainty into the model, <span> \( x_{t+1} = z_g f(x_t) \) </span>, where `z_g` is lognormal with logmean 0 and logsd of `0.15`.  








As before, we simulate 100 replicates using the same random number sequence, now under this case where the noise in growth is intrinsic and is being accounted for by the management.  




### Growth uncertainty & bias  





This time we consider the same optimization under uncertainty as before, but the simulations introduce bias through a random estimate of the growth rate parameter A, drawn from a normal with mean equal to the true value `1.5` and variance `0.2`.   Since A is a constant multiplier in the growth dynamics, this is equivalent to a random estimate of mean of the noise process `z_g`.  Our estimate of the parameter is drawn from the distribution and then held fixed for that replicate.  Each replicate draws its own value, so average parameter estimate across the replicates should be close to the true value.  _Isn't this equivalent to the standard parameter uncertainty?_




For the record, the biases by index number are:



```
  [1] 1.766 1.380 1.511 1.394 1.484 1.532 1.480 1.429 1.406 1.507 1.460
 [12] 1.391 1.575 1.581 1.534 1.280 1.580 1.585 1.443 1.611 1.753 1.566
 [23] 1.653 1.757 1.832 1.404 1.418 1.434 1.720 1.616 1.692 1.560 1.442
 [34] 1.209 1.159 1.902 1.501 1.612 1.639 1.789 1.627 1.505 1.649 1.569
 [45] 1.374 1.692 1.533 1.493 1.316 1.730 1.422 1.480 1.494 1.489 1.611
 [56] 1.680 1.607 1.206 1.018 1.427 1.464 1.470 1.221 1.364 1.136 2.028
 [67] 1.609 1.499 1.122 1.447 1.664 1.656 1.591 1.782 1.138 1.280 1.496
 [78] 1.530 1.394 1.740 1.091 1.297 1.821 1.678 1.196 1.665 1.573 1.234
 [89] 1.572 1.679 1.301 1.604 1.253 1.606 1.303 1.425 1.292 1.640 1.471
[100] 1.600
```




## Known bias

The correct baseline comparison against the above scenario, Mike points out to me, still includes the bias, but solves for the optimal solution of that random draw.  This is much slower since it requires solving a new optimum each time, but parallelizes well.   







## Summarize and plot the results                                                   





### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  

![plot of chunk onerep](http://farm8.staticflickr.com/7241/7223268132_a05f71a08f_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7227/7223268766_182021509e_o.png) 



![The profits made in each time interval of a single replicate, by scenario](http://farm8.staticflickr.com/7100/7223269436_3e762aa435_o.png) 



![the distribution of profits by scenario](http://farm9.staticflickr.com/8146/7223269932_f557f2b5e4_o.png) 


Summary statistics 



```
     uncertainty    V1
[1,]       known 33.08
[2,]      Growth 34.34
[3,]  RandomBias 35.47
[4,]   KnownBias 35.47
```



```
     uncertainty     V1
[1,]       known  0.000
[2,]      Growth  3.963
[3,]  RandomBias 15.403
[4,]   KnownBias 15.403
```








Direct comparisons: 



```r
setkey(dt, uncertainty)
cost_of_bias = dt[uncertainty=="KnownBias", sum(profit), by=reps]$V1 - dt[uncertainty=="RandomBias", sum(profit), by=reps]$V1
cost_of_bias
```



```
  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [36] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [71] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```



```r
mean(cost_of_bias)
```



```
[1] 0
```



```r
sd(cost_of_bias)
```



```
[1] 0
```




# References

Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005). "Fishery
management under multiple uncertainty." _Journal of Environmental
Economics and Management_, *50*. ISSN 00950696, <URL:
http://dx.doi.org/10.1016/j.jeem.2004.11.005>.


