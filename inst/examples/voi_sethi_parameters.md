






# Calculating the value of information

 Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).
 Compute the optimal solution under different forms of uncertainty.   










Chose the state equation / population dynamics function



```r
f <- function(x,h,p){
	S = x - h
	p[1] * S * (1 - S/p[2]) + S
}
```




With parameters `r` = `1` and `K` = `100`.



```r
pars <- c(r, K)
```




We consider a profits from fishing to be a function of harvest `h` and stock size `x`,  \\( \Pi(x,h) = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\),



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




We seek a harvest policy which maximizes the discounted profit from the fishery using a stochastic dynamic programming approach over a discrete grid of stock sizes from `0` to `150` on a grid of `100` points, and over an identical discrete grid of possible harvest values.  




```r
x_grid <- seq(xmin, xmax, length = grid_n)  
h_grid <- x_grid  
```







```r
delta <- 0.05
xT <- 0
OptTime <- 25
```




We will determine the optimal solution over a `25` time step window with boundary condition for stock at `0` and discounting rate of `0.05`.  

# Scenarios: 

We use Monte Carlo integration over the noise processes to determine the transition matrix.  



```r
require(snowfall)
sfInit(parallel=TRUE, cpu=16)
```

```
R Version:  R version 2.14.1 (2011-12-22) 

```







```r
scenario <- function(policy_g, policy_m, policy_i){ 

  z_g <- function() 1+(2*runif(1, 0,  1)-1) * policy_g
  z_m <- function() 1+(2*runif(1, 0,  1)-1) * policy_m
  z_i <- function() 1+(2*runif(1, 0,  1)-1) * policy_i

  SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=2e4)
  opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                       profit, delta, reward=0)
}
```






```r
det <- scenario(0.1, 0.0, 0.0)
```

```
Library ggplot2 loaded.
```

```r
g   <- scenario(0.5, 0.0, 0.0)
```

```
Library ggplot2 loaded.
```

```r
m   <- scenario(0.0, 0.5, 0.0)
```

```
Library ggplot2 loaded.
```

```r
i   <- scenario(0.0, 0.0, 0.5)
```

```
Library ggplot2 loaded.
```

```r
gm  <- scenario(0.5, 0.5, 0.0)
```

```
Library ggplot2 loaded.
```

```r
gi  <- scenario(0.5, 0.0, 0.5)
```

```
Library ggplot2 loaded.
```

```r
im  <- scenario(0.0, 0.5, 0.5)
```

```
Library ggplot2 loaded.
```

```r
gmi <- scenario(0.5, 0.5, 0.5)
```

```
Library ggplot2 loaded.
```





### plots



```r
require(reshape2)
policy <- melt(data.frame(stock = x_grid, det = det$D[,1], g = g$D[,1], m = m$D[,1], i= m$D[,i], gm = gm$D[,1], gi = gi$D[,1], gmi = gmi$D[,1]), id = "stock")
```

```
Error: invalid subscript type 'list'```

```r
ggplot(policy) + geom_point(aes(stock, stock-x_grid[value], color=variable)) 
```

```
Error: object 'policy' not found```

```r
	geom_smooth(aes(stock, stock-x_grid[value], color=variable)) + ylab("escapement") 
```

```
Error: non-numeric argument to binary operator```

```r

ggplot(policy) + geom_point(aes(stock, stock-x_grid[value], color=variable)) 
```

```
Error: object 'policy' not found```

```r
	geom_smooth(aes(stock, x_grid[value], color=variable)) + ylab("harvest") 
```

```
Error: non-numeric argument to binary operator```

```r


value <- melt(data.frame(stock = x_grid, det=det$V, g = g$V, m = m$V, gm = gm$V, gi = gi$V, gmi = gmi$V), id = "stock")
ggplot(value) + geom_point(aes(stock, value, color=variable)) +
	ggplot(value) + geom_smooth(aes(stock, value, color=variable)) + ylab("Net Present Value")
```

```
Error: non-numeric argument to binary operator```




## Simulations



```r
simulatereps <- function(opt, true_g, true_m, true_i){
  z_g <- function() rlnorm(1,  0, true_g) 
  z_m <- function() rlnorm(1,  0, true_m) 
  z_i <- function() rlnorm(1,  0, true_i) 

  sims <- lapply(1:100, function(i){
    ForwardSimulate(f, pars, x_grid, h_grid, x0 = K, opt$D, z_g, z_m, z_i, profit)
  })
  
#  list(sims = sims, opt = opt, true_stochasticity = c(true_g, true_m, true_i))
  sims
}
```





All cases



```r
policyfn <- list(det=det, g=g, m=m, i=i, gm=gm, gi=gi, gmi=gmi)
noise <- list(s0=c(0,0,0), sg=c(.2, 0, 0), sm=c(0, .2, 0), si=c(0,0,.2), sgm=c(.2, .2, 0), sgi=c(.2, 0, .2), sgmi=c(.2, .2, .2)) 
allcases <- lapply(policyfn, function(policyfn_i){
  lapply(noise, function(noise_i){
    simulatereps(policyfn_i, noise_i[1], noise_i[2], noise_i[3])
  })
})
```






```r
sims <- unlist(allcases, recursive=FALSE)
dat <- melt(sims, id=names(sims[[1]][[1]]))  
dt <- data.table(dat)
setnames(dt, c("L2", "L1"), c("reps", "uncertainty")) # names are nice
```





### Plots 




```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_line(aes(time, harvest), col="darkgreen") + 
  facet_wrap(~uncertainty) 
```

![plot of chunk onerep](http://farm8.staticflickr.com/7240/7398748034_a50fe77636_o.png) 


This plot summarizes the stock dynamics by visualizing the replicates.



```r
p1 <- ggplot(dt) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm6.staticflickr.com/5276/7398751606_bb702f69d6_o.png) 





```r
profits <-dt[ , sum(profit), by=c("reps", "uncertainty")] 
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7228/7398753726_16845bdbc8_o.png) 


Summary statistics 



```r
means <- profits[, mean(V1), by=uncertainty]
sds <- profits[, sd(V1), by=uncertainty]
```






```r
require(xtable)
uncertainties <- c("deterministic", "growth", "measure", "implement", "growth+measure", "growth+implement", "measure+implement", "all")
print(xtable(matrix(means$V1, nrow=7, dimnames=list(uncertainties, uncertainties))), type="html")
```

```
Error: length of 'dimnames' [1] not equal to array extent```

```r
print(xtable(matrix(sds$V1, nrow=7, dimnames=list(uncertainties, uncertainties))), type="html")
```

```
Error: length of 'dimnames' [1] not equal to array extent```






# References

