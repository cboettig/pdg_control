






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





This time we consider the same optimization under uncertainty as before, but the simulations introduce bias through a random estimate of the growth rate parameter A, drawn from a normal with mean equal to the true value `1.5` and variance `0.1`.   Since A is a constant multiplier in the growth dynamics, this is equivalent to a random estimate of mean of the noise process `z_g`.  Our estimate of the parameter is drawn from the distribution and then held fixed for that replicate.  Each replicate draws its own value, so average parameter estimate across the replicates should be close to the true value.  _Isn't this equivalent to the standard parameter uncertainty?_





## Known bias

The correct baseline comparison against the above scenario, Mike points out to me, still includes the bias, but solves for the optimal solution of that random draw.  This is much slower since it requires solving a new optimum each time, but parallelizes well.   




```
R Version:  R version 2.14.1 (2011-12-22) 

```



```
Library pdgControl loaded.
```






## Summarize and plot the results                                                   





### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  

![plot of chunk onerep](http://farm9.staticflickr.com/8016/7223092804_9ef9492b06_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7228/7223093476_421cb98440_o.png) 



![The profits made in each time interval of a single replicate, by scenario](http://farm6.staticflickr.com/5470/7223094000_419a639c96_o.png) 



![the distribution of profits by scenario](http://farm8.staticflickr.com/7241/7223094524_2463b1a39a_o.png) 


Summary statistics 



```
     uncertainty    V1
[1,]       known 33.08
[2,]      Growth 34.34
[3,]  RandomBias 34.66
[4,]   KnownBias 34.86
```



```
     uncertainty    V1
[1,]       known 0.000
[2,]      Growth 3.963
[3,]  RandomBias 9.497
[4,]   KnownBias 9.648
```








Direct comparisons: 



```r
setkey(dt, uncertainty)
cost_of_bias = dt[uncertainty=="KnownBias", sum(profit), by=reps]$V1 - dt[uncertainty=="RandomBias", sum(profit), by=reps]$V1
cost_of_bias
```



```
  [1]  11.6827   2.5336  -4.3163   3.8347  -0.8089  -2.0183   1.2863
  [8]  -4.1427  -4.3798   4.4432   3.4784 -11.2389  -1.7576   8.8768
 [15]  -0.2049   3.5902   0.1604  -5.5709  -6.2378  15.6252   5.0021
 [22]   1.1942   3.9689   1.2102   6.0537  -1.8033   4.9185   1.4861
 [29]   5.4467   0.1510   9.4859  -3.8346  -6.8022  -5.0448  -0.8049
 [36]  -3.0257  -0.7626  -5.3635   3.6319   5.7558   4.0046  -1.8527
 [43]  -8.9348   0.6039  -1.6645   6.9408   1.2359  -6.3933   6.6986
 [50] -10.2909   0.4047   2.3933   0.7641  -2.7136   0.8972   5.6496
 [57]   8.6733   4.0855   0.7901  -3.0276   6.0774  -4.7135 -12.3060
 [64]   0.6089  -2.6245  -7.3993  -2.8139  11.9580   9.6845  -6.0542
 [71]   9.9349  -7.8665   1.2105 -10.3859  -3.0245  -2.6540   8.6623
 [78] -14.0362  -0.7954  -0.5253  -5.3299  -0.4057   7.8684  10.2892
 [85]  -3.6322  -3.8366  -3.8148   1.4110  -5.7724   9.6865   2.6215
 [92]  -4.9402  -1.5327  -7.9145  -1.0575   4.0332  -2.1267  -4.6049
 [99]   9.0792   2.7995
```



```r
mean(cost_of_bias)
```



```
[1] 0.1972
```



```r
sd(cost_of_bias)
```



```
[1] 5.849
```




# References

Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005). "Fishery
management under multiple uncertainty." _Journal of Environmental
Economics and Management_, *50*. ISSN 00950696, <URL:
http://dx.doi.org/10.1016/j.jeem.2004.11.005>.


