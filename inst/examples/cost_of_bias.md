






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

We calculate the stochastic transition matrix for the probability of going from any state \\(x_t \\) to any other state \\(x_{t+1}\\) the following year, for each possible choice of harvest \\( h_t \\).  This provides a look-up table for the dynamic programming calculations. In the Sethi case, computing the distribution over multiple sources of noise is actually quite difficult.  Simulation turns out to be more efficient than numerically integrating over each distribution.  We use this matrix to compute the optimum strategy for all possible states of the world by dynamic programming, and then simulate replicates while applying this rule.   


### No Uncertainty 

The first scenario considers the completely deterministic case.  










We simulate 100 replicates of this system.  We will used a fixed seed so that we can compare these replicates to simulations under different conditions.  (Of course the seed is irrelevant at this stage since this is actually deterministic).  



### Growth uncertainty 




The next scenario introduces growth uncertainty into the model, <span> \( x_{t+1} = z_g f(x_t) \) </span>, where `z_g` is lognormal with logmean 0 and logsd of `0.15`.  








As before, we simulate 100 replicates using the same random number sequence, now under this case where the noise in growth is intrinsic and is being accounted for by the management.  




### Growth uncertainty & bias  


`` {r }
est_pars <- pars
par_var <- .1
````

This time we consider the same optimization under uncertainty as before, but the simulations introduce bias through a random estimate of the growth rate parameter A, drawn from a normal with mean equal to the true value `1.5` and variance `0.1`.   Since A is a constant multiplier in the growth dynamics, this is equivalent to a random estimate of mean of the noise process `z_g`.  Our estimate of the parameter is drawn from the distribution and then held fixed for that replicate.  Each replicate draws its own value, so average parameter estimate across the replicates should be close to the true value.  _Isn't this equivalent to the standard parameter uncertainty?_





## Unbiased error in growth estimate




We will compare this result to the situation of underestimating the uncertainty in the growth rate, but knowning the parameter exactly.  We use the deterministic optimum solution under a reality that has log-normal growth noise equal to the sum of the variances in the previous example, \( \sigma_g = \)  `0.1803`.  







## Summarize and plot the results                                                   





### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  

![plot of chunk onerep](http://farm9.staticflickr.com/8008/7212181630_8f1805db9b_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm6.staticflickr.com/5466/7212182138_d38efebe87_o.png) 



![The profits made in each time interval of a single replicate, by scenario](http://farm8.staticflickr.com/7092/7212182474_e90a5663a1_o.png) 



![the distribution of profits by scenario](http://farm8.staticflickr.com/7216/7212182814_af80a6f05e_o.png) 


Summary statistics tell the final story:



```
       uncertainty    V1
[1,]         known 33.06
[2,]        Growth 34.34
[3,]          Bias 35.38
[4,] Underestimate 34.92
```



```
       uncertainty    V1
[1,]         known 0.000
[2,]        Growth 3.963
[3,]          Bias 8.680
[4,] Underestimate 4.912
```






# References

Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
"Fishery management under multiple uncertainty." _Journal of
Environmental Economics and Management_, *50*. ISSN 00950696,
<URL: http://dx.doi.org/10.1016/j.jeem.2004.11.005>.


