# file plots.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-16
# creates extra plots accompanying Reed.R

# for stat plots
require(ggplot2)
require(Hmisc)
require(pdgControl)

## FIXME Once standardized, all these plots should become package fns

## Reshape and summarize data ###
dat <- melt(sims, id="time") # reshapes the data matrix to "long" form
# some stats on the replicates, (stat_sumary can do this instead)
m <- cast(dat, time ~ variable, mean) # mean population
err <- cast(dat, time ~ variable, sd) # sd population


## Show dynamics of a single replicate 
ex <- sample(1:100,1) # a random replicate
p0 <- ggplot(dat) +
      geom_line(aes(time, value, color = variable), 
                data = subset(dat, L1 == ex)) + 
      geom_abline(intercept = opt$S, slope = 0, col = "darkred") + # show Reed's S: optimal escapement 
      geom_abline(intercept = xT, slope = 0, lty = 2) #+ # show Allee threshold
#print(p0)

p1 <- plot_replicates(sims)

## Update the main p1 plot to 
## make the crashed trajectories stand out?
#p1 <- p1 + geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "harvest" & 
#          (L1 %in% optimal_crashed$L1)),  col = "darkgreen", alpha = 0.5) +
#     geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "fishstock" &
#          (L1 %in% optimal_crashed$L1)),  col = "darkblue", alpha = 0.5) 


#### p2 Plots the unharvested dynamics ###
crashed = subset(dat, variable =="unharvested" & time == OptTime - 1 & value < xT)
p2 <- ggplot(data = subset(dat, variable == "unharvested"),
             aes(time, value, group = L1)) + geom_line(alpha = 0.2) + 
  # shows the mean +/- mult * sd , requires Hmisc
  stat_summary(mapping = aes(group = 1), fun.data = mean_sdl, 
               geom = "smooth", mult = 1) +
  opts(title=sprintf("Unfished dynamics, %d populations crash", dim(crashed)[1]))



## Profits plot #######
p3 <- ggplot(subset(dat, variable == "harvest"), 
             aes(time, profit(value, K), group = L1)) + 
  geom_line(alpha = .2) + 
  labs(x = "Time (yrs)", y = "Profit") + 
  stat_summary(fun.data = mean_sdl, 
    geom = "smooth", mapping = aes(group = 1),
    lwd = 1, col = "darkred", mult = 1)

# fun.data should be used for functions that give mean+sd back
# otherwise, use fun.y = mean, fun.ymin = 

# profits for each rep, by timestep
require(plyr)
cash <- ddply(subset(dat, variable == "harvest"), "L1", 
  function(df) profit(dat$variable, K))
require(data.table)
DT <- data.table(dat)
DT[, profit

# add mean line by hand
# p3 <- p3+geom_line(aes(time,rowMeans(cash))) # average profit made as function of time


p4 <- qplot(colSums(cash), xlab="Total Profit", ylab=NULL) + # histogram of total profit made
  geom_vline(xintercept=mean(colSums(cash))) + # expected total profit
  opts(plot.margin=unit(rep(0,4), "lines")) + theme_gray(9) #appearances
subvp <- viewport(width=.3, height=.3, x=.8, y=.8)

#png("profit.png")
#print(p3)
#print(p4, vp=subvp)
#dev.off()
#ggsave("samplerun.png", plot=p0)
#ggsave("fished.png", plot=p1)
#ggsave("unfished.png", plot=p2)

#require(socialR)
#upload("profit.png samplerun.png fished.png unfished.png", script="Reed.R", tag="PDG_Control")
