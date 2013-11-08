

```r
rm(list = ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```






Define the parameters of our cost function. We will 
consider a fixed price and fixed $c_0$ cost per unit effort.  
An additional "policy cost" is introduced through the $c_2$
coefficient, which can take a range of values. We will loop
over this range for each of the functional forms of the policy
cost, in order to choose coefficients $c_2$ in each case
that are comparable.  



```r
price = 10
c0 = 30
profit <- profit_harvest(price = price, c0 = c0, c1 = 0)
```



```r
c2 <- exp(seq(0, log(41), length.out = 40)) - 1
c2 <- seq(0, 40, length.out = 100)
```



```r
reduction <- 0.25
```


We will compare when transaction costs have reduced the value to `(1-reduction) * adjustment_cost_free_value`




```r
opts_knit$set(upload.fun = socialR::flickr.url)
opts_chunk$set(dev.args = list(bg = "transparent"), tidy = FALSE, comment = NA, 
    message = FALSE, warning = FALSE)
opts_chunk$set(echo = FALSE)
theme_set(theme_bw())
```




This block defines the various parameters and nuisance parameters we
need to specify our optimal control problem.






Given these parameters, we can determine the optimal solution under the
classic assumption of "adjustment-free" costs, where there is no penalty
for adjusting the value of the control (harvest) at each decision step
(year).






Now we introduce the functional forms for each of the adjustment costs,






## Apples to Apples levels

Before we can start comparing solutions induced under each of these
costs, we need to scale them appropriately (e.g. the coefficient $c_2$
has different units and a different magnitude of effect when multiplied
by quadratic differences $(h_{t} - h_{t-1})^2$ then when multiplied by
linear differences $h_{t} - h_{t-1}$.


To do so we loop over a grid of $c_2$ values and determine the net
present value under each functional form of adjustment costs for each
$c_2$ in our grid.



```
R Version:  R version 3.0.2 (2013-09-25) 
```

```
Library pdgControl loaded.
```



### Loop over penalty functions and magnitudes






Note that `optim_policy` has been updated to return the equilibrium value
of profits from fish harvests before the adjustment costs have been paid,
`penalty_free_V`, as well as the total value `V`. Initially we based the 
comparison on matching `penalty_free_V`, now the economists recommend we
simply use the total net present value resulting, `V`.  



The value matrices returned in either case contain the values for all possible states,
we simply evaluate it at the carrying capacity (which is our initial
condition.)  The index in `x_grid` that corresponds to the carrying
capacity (initial condition) `i` indicates this.



Extract the policy cost 





Tidy up the data and plot the net present value relative to that achieved when managed without a penalty.


```
[1] 190.6
```

![plot of chunk npv-plot](http://farm8.staticflickr.com/7300/10736952854_044a56201f_o.png) 


Alternative version of plot, showing a ratio

![plot of chunk apples_plot](http://farm3.staticflickr.com/2816/10736948206_94417d773b_o.png) 



Find the value of `c2` that brings each penalty closest to 25% of the cost-free adjustment value:


```
   L2    L1 fixed 
0.404 1.616 9.293 
```







How close has each of these actually gotten to 25% reduction to NPV0 (which is 142.9516)?




<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Thu Nov  7 21:16:04 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> model </TH> <TH> value.realized </TH> <TH> percent.of.npv0 </TH> <TH> percent.error </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> L2 </TD> <TD align="right"> 152.42 </TD> <TD align="right"> 79.97 </TD> <TD align="right"> 6.62 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> L1 </TD> <TD align="right"> 142.37 </TD> <TD align="right"> 74.69 </TD> <TD align="right"> -0.41 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> fixed </TD> <TD align="right"> 143.10 </TD> <TD align="right"> 75.08 </TD> <TD align="right"> 0.11 </TD> </TR>
   </TABLE>



Solve the policy cost for the specified penalty function

















```
  model value.realized percent.of.npv0 percent.error
1    L2          152.4           79.97        6.6244
2    L1          142.4           74.69       -0.4091
3 fixed          143.1           75.08        0.1052
```

```
[1] -1.265e+07  7.510e-01 -1.265e+07
```









### additional plots based on simulations under the correct model





```
          L1     L2    fixed
var 7.116242 5.3010  8.91387
a   0.007221 0.2078 -0.04762
```




# Plots 


![plot of chunk p1](http://farm4.staticflickr.com/3818/10737166943_ec84cf11ab_o.png) 


![plot of chunk Figure3](http://farm4.staticflickr.com/3748/10736949806_9cab36457c_o.png) 


# Figure 4





















![plot of chunk Figure4](http://farm8.staticflickr.com/7395/10736954766_47e1d02091_o.png) ![plot of chunk Figure4](http://farm4.staticflickr.com/3690/10737172953_d5d31cd62f_o.png) 



![plot of chunk Figure4S](http://farm4.staticflickr.com/3706/10736955576_00ca28919d_o.png) ![plot of chunk Figure4S](http://farm4.staticflickr.com/3733/10737173753_0b07810c44_o.png) 










