



```r
require(pdgControl)
```

```
## Loading required package: pdgControl
```

```r
profit <- profit_harvest(price = 10, c0 = 0, c1 = 0)
c2 <- exp(seq(0, log(21), length.out = 20)) - 1
```



```r
knit("policycosts.Rmd", "policycosts_c0_0.md")
```

```
## 
## 
## processing file: policycosts.Rmd
```

```
## output file:
## /home/cboettig/Computing/pdgControl/inst/examples/policycosts/policycosts_c0_0.md
```


## Results


```r
ggplot(dat, aes(c2, (npv0 - value)/npv0, col = variable)) + geom_point() + geom_line()
```

```
## Error: object 'npv0' not found
```



```r
apples
```

```
##            L2            L1 free_decrease         fixed free_increase 
##       11.9852        7.0293        0.0000       20.0000        0.0000 
##          quad 
##        0.8983
```




```r
p0
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-12-03-3ee87d3584-unnamed-chunk-4.png) 



```r
p1
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-03-3ee87d3584-unnamed-chunk-51.png) 

```r
p2
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-03-3ee87d3584-unnamed-chunk-52.png) 



