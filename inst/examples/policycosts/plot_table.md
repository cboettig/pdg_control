

```r
opts_knit$set(upload.fun = socialR::flickr.url)
```



```r
require(ggplot2)
df <- read.csv("output.csv")
df <- df[2:8]

require(reshape2)
who <- c("penalty_fn", "ignore_fraction", "assume_fraction", "sigma_g", "reduction")
df2 <- melt(df, id = who)
who <- c("penalty_fn", "sigma_g", "reduction", "variable", "value")
df3 <- melt(df2, id = who)
names(df3) <- c("penalty_fn", "sigma_g", "reduction", "variable", "value", "fraction_variable", 
    "fraction")

df3$reduction <- as.factor(df3$reduction)
# ggplot(df3) + geom_line(aes(sigma_g, value, col=penalty_fn, lty=variable))
# + facet_wrap(~reduction)
ggplot(df3) + geom_line(aes(sigma_g, value, col = reduction, lty = variable)) + 
    facet_wrap(~penalty_fn)
```

![plot of chunk unnamed-chunk-2](http://farm6.staticflickr.com/5495/10693194473_551b21647d_o.png) 




```r
print(xtable::xtable(df3), type = "html")
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Tue Nov  5 12:12:14 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> penalty_fn </TH> <TH> sigma_g </TH> <TH> reduction </TH> <TH> variable </TH> <TH> value </TH> <TH> fraction_variable </TH> <TH> fraction </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 14536.09 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 11020.10 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.76 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 13584.97 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 9273.66 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.09 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 11020.10 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.76 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 8332.45 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.84 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -641.78 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.21 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -272586.11 </TD> <TD> ignore_fraction </TD> <TD align="right"> -24.90 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -1213.01 </TD> <TD> ignore_fraction </TD> <TD align="right"> -0.46 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 18281.97 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.01 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 12462.51 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.72 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 16691.94 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.03 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 12788.95 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.96 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 12462.51 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.72 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 12399.01 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.88 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 9099.48 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.73 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -6337.76 </TD> <TD> ignore_fraction </TD> <TD align="right"> -0.39 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 7399.01 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.62 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 30227.43 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.07 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 33472.18 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 28926.73 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.09 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 24157.93 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.10 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 19091.31 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.93 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 23573.19 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.02 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 21086.12 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.12 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 19091.31 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.93 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 17209.56 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.85 </TD> </TR>
  <TR> <TD align="right"> 28 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 16857.61 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 29 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 15538.44 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.76 </TD> </TR>
  <TR> <TD align="right"> 30 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 17582.78 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 31 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 17561.81 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.09 </TD> </TR>
  <TR> <TD align="right"> 32 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 15538.44 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.76 </TD> </TR>
  <TR> <TD align="right"> 33 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 17401.82 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.84 </TD> </TR>
  <TR> <TD align="right"> 34 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 15451.25 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.21 </TD> </TR>
  <TR> <TD align="right"> 35 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 11567.12 </TD> <TD> ignore_fraction </TD> <TD align="right"> -24.90 </TD> </TR>
  <TR> <TD align="right"> 36 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 15519.95 </TD> <TD> ignore_fraction </TD> <TD align="right"> -0.46 </TD> </TR>
  <TR> <TD align="right"> 37 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 20973.56 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.01 </TD> </TR>
  <TR> <TD align="right"> 38 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 19392.85 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.72 </TD> </TR>
  <TR> <TD align="right"> 39 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 20143.80 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.03 </TD> </TR>
  <TR> <TD align="right"> 40 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 20663.21 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.96 </TD> </TR>
  <TR> <TD align="right"> 41 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 19392.85 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.72 </TD> </TR>
  <TR> <TD align="right"> 42 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 21289.21 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.88 </TD> </TR>
  <TR> <TD align="right"> 43 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 19999.07 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.73 </TD> </TR>
  <TR> <TD align="right"> 44 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 17765.33 </TD> <TD> ignore_fraction </TD> <TD align="right"> -0.39 </TD> </TR>
  <TR> <TD align="right"> 45 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 21969.84 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.62 </TD> </TR>
  <TR> <TD align="right"> 46 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 31377.07 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.07 </TD> </TR>
  <TR> <TD align="right"> 47 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 33472.18 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 48 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 30500.85 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.09 </TD> </TR>
  <TR> <TD align="right"> 49 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 31767.92 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.10 </TD> </TR>
  <TR> <TD align="right"> 50 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 29225.46 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.93 </TD> </TR>
  <TR> <TD align="right"> 51 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 30069.80 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.02 </TD> </TR>
  <TR> <TD align="right"> 52 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 32035.93 </TD> <TD> ignore_fraction </TD> <TD align="right"> 1.12 </TD> </TR>
  <TR> <TD align="right"> 53 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 29225.46 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.93 </TD> </TR>
  <TR> <TD align="right"> 54 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 32946.24 </TD> <TD> ignore_fraction </TD> <TD align="right"> 0.85 </TD> </TR>
  <TR> <TD align="right"> 55 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 14536.09 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 56 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 11020.10 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 57 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 13584.97 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.05 </TD> </TR>
  <TR> <TD align="right"> 58 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 9273.66 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.04 </TD> </TR>
  <TR> <TD align="right"> 59 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 11020.10 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 60 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 8332.45 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.03 </TD> </TR>
  <TR> <TD align="right"> 61 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -641.78 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 62 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -272586.11 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.69 </TD> </TR>
  <TR> <TD align="right"> 63 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -1213.01 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 64 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 18281.97 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.99 </TD> </TR>
  <TR> <TD align="right"> 65 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 12462.51 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> 66 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 16691.94 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.95 </TD> </TR>
  <TR> <TD align="right"> 67 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 12788.95 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.97 </TD> </TR>
  <TR> <TD align="right"> 68 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 12462.51 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> 69 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 12399.01 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 70 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 9099.48 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.94 </TD> </TR>
  <TR> <TD align="right"> 71 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> -6337.76 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.83 </TD> </TR>
  <TR> <TD align="right"> 72 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 7399.01 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.03 </TD> </TR>
  <TR> <TD align="right"> 73 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 30227.43 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.94 </TD> </TR>
  <TR> <TD align="right"> 74 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 33472.18 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 75 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> ignore_cost </TD> <TD align="right"> 28926.73 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> 76 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 24157.93 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.95 </TD> </TR>
  <TR> <TD align="right"> 77 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 19091.31 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.87 </TD> </TR>
  <TR> <TD align="right"> 78 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> ignore_cost </TD> <TD align="right"> 23573.19 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.90 </TD> </TR>
  <TR> <TD align="right"> 79 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 21086.12 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.96 </TD> </TR>
  <TR> <TD align="right"> 80 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 19091.31 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.87 </TD> </TR>
  <TR> <TD align="right"> 81 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> ignore_cost </TD> <TD align="right"> 17209.56 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.98 </TD> </TR>
  <TR> <TD align="right"> 82 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 16857.61 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 83 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 15538.44 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 84 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 17582.78 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.05 </TD> </TR>
  <TR> <TD align="right"> 85 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 17561.81 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.04 </TD> </TR>
  <TR> <TD align="right"> 86 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 15538.44 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 87 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 17401.82 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.03 </TD> </TR>
  <TR> <TD align="right"> 88 </TD> <TD> L1 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 15451.25 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 89 </TD> <TD> L2 </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 11567.12 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.69 </TD> </TR>
  <TR> <TD align="right"> 90 </TD> <TD> fixed </TD> <TD align="right"> 0.05 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 15519.95 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.92 </TD> </TR>
  <TR> <TD align="right"> 91 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 20973.56 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.99 </TD> </TR>
  <TR> <TD align="right"> 92 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 19392.85 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> 93 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 20143.80 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.95 </TD> </TR>
  <TR> <TD align="right"> 94 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 20663.21 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.97 </TD> </TR>
  <TR> <TD align="right"> 95 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 19392.85 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> 96 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 21289.21 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 97 </TD> <TD> L1 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 19999.07 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.94 </TD> </TR>
  <TR> <TD align="right"> 98 </TD> <TD> L2 </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 17765.33 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.83 </TD> </TR>
  <TR> <TD align="right"> 99 </TD> <TD> fixed </TD> <TD align="right"> 0.20 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 21969.84 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.03 </TD> </TR>
  <TR> <TD align="right"> 100 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 31377.07 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.94 </TD> </TR>
  <TR> <TD align="right"> 101 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 33472.18 </TD> <TD> assume_fraction </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 102 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.1 </TD> <TD> assume_cost </TD> <TD align="right"> 30500.85 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.91 </TD> </TR>
  <TR> <TD align="right"> 103 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 31767.92 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.95 </TD> </TR>
  <TR> <TD align="right"> 104 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 29225.46 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.87 </TD> </TR>
  <TR> <TD align="right"> 105 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.2 </TD> <TD> assume_cost </TD> <TD align="right"> 30069.80 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.90 </TD> </TR>
  <TR> <TD align="right"> 106 </TD> <TD> L1 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 32035.93 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.96 </TD> </TR>
  <TR> <TD align="right"> 107 </TD> <TD> L2 </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 29225.46 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.87 </TD> </TR>
  <TR> <TD align="right"> 108 </TD> <TD> fixed </TD> <TD align="right"> 0.50 </TD> <TD> 0.3 </TD> <TD> assume_cost </TD> <TD align="right"> 32946.24 </TD> <TD> assume_fraction </TD> <TD align="right"> 0.98 </TD> </TR>
   </TABLE>

