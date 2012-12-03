



```r
profit <- profit_harvest(price = 10, c0 = 30, c1 = 0)
```

```
## Error: could not find function "profit_harvest"
```

```r
c2 <- seq(0, 20, length.out = 40)
```



```r
knit("policycosts.Rmd", "policycosts_c0_30.md")
```

```
## 
## 
## processing file: policycosts.Rmd
```

```
## output file:
## /home/cboettig/Computing/pdg_control/inst/examples/policycosts/policycosts_c0_30.md
```


## Results


```r
ggplot(dat, aes(c2, (npv0 - value)/npv0, col = variable)) + geom_point() + geom_line()
```

```
## Error: object 'value' not found
```



```r
as.table(apples)
```

```
##            L2            L1 free_decrease         fixed free_increase 
##         9.231         5.641         0.000        15.385         0.000 
##          quad 
##         4.103
```




```r
p0
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-12-03-2d28efb6b5-unnamed-chunk-4.png) 



```r
p1
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-03-2d28efb6b5-unnamed-chunk-51.png) 

```r
p2
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-03-2d28efb6b5-unnamed-chunk-52.png) 



