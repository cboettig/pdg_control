# file plots.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-16
# creates extra plots accompanying Reed.R


## Example plot the results of a single run, against unharvested version  
out <- ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D,
                       z_g, z_m, z_i, interval=1)
dat <- melt(out, id="time")
p0 <- ggplot(dat, aes(time, value, color=variable)) + geom_line() +  
  geom_abline(intercept=opt$S, slope=0, col="black") + # Reed's S,
  geom_abline(intercept = xT, slope=0, col="darkred", lty=3) # unfished Allee 
p0 <- p0 + geom_abline(intercept=e_star, slope=0, col="green", lty = 2) # tippt


## Update the main p1 plot to 
## make the crashed trajectories stand out?
#p1 <- p1 + geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "harvest" & 
#          (L1 %in% optimal_crashed$L1)),  col = "darkgreen", alpha = 0.5) +
#     geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "fishstock" &
#          (L1 %in% optimal_crashed$L1)),  col = "darkblue", alpha = 0.5) 


crashed = subset(dat, variable =="unharvested" & time == OptTime-1 & value < xT)
p2 <- ggplot(dat) + 
  geom_line(aes(time, value, group = L1), 
            data = subset(dat, variable == "unharvested"), alpha=.2) + 
  geom_line(aes(time, cast(dat, time ~ variable, mean)$unharvested)) 

p2 <- p2 + opts(title=sprintf("Unfished dynamics, %d populations crash",
  dim(crashed)[1]))


## Profits plot
p3 <- ggplot(subset(dat, variable == "harvest"), aes(profit(x_grid,value))) + geom_histogram()

# calculate end profits:
d <- subset(dat, variable == "harvest") 
income <- profit(x_grid,d$value)

print(p3)


#ggsave("samplerun.png", plot=p0)
#ggsave("fished.png", plot=p1)
#ggsave("unfished.png", plot=p2)


# check out the deterministic dynamics
x <- numeric(OptTime)
x[1] <- K
for(t in 1:(OptTime-1))
  x[t+1] <- f(x[t],e_star+.01,pars)


