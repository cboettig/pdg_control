# file plots.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-16
# creates extra plots accompanying Reed.R

# for stat plots
require(Hmisc)

## Show dynamics of a single replicate 
ex <- sample(1:100,1) # a random replicate
p0 <- ggplot(dat) +
      geom_line(aes(time, value, color=variable), 
                data = subset(dat, L1 == ex)) + 
      geom_abline(intercept=opt$S, slope = 0, col="darkred") + # show Reed's S: optimal escapement 
      geom_abline(intercept=xT, slope = 0, lty=2) #+ # show Allee threshold
print(p0)



## Show the ensemble fishstock and harvest dynamics 
p1 <- ggplot(dat) +
      # Replicate harvested dynamics
      geom_line(aes(time, value, group = L1), data = 
                subset(dat, variable == "fishstock"), alpha = 0.2) + 
      ## Mean & SD for population
      geom_ribbon(aes(x = time, ymin = m$fishstock - err$fishstock,
                      ymax = m$fishstock + err$fishstock),
                  fill = "darkblue", alpha = 0.4)  +
      geom_line(aes(time, m$fishstock), col = "lightblue")  +
      geom_abline(intercept=opt$S, slope = 0) + # show Reed's S: optimal escapement 
      geom_abline(intercept=xT, slope = 0, lty=2) + # show Allee threshold
      ## And the same for harvest-effort levels
      geom_line(aes(time, value, group = L1), 
            data = subset(dat, variable == "harvest"),  alpha=.2, col = "darkgreen") +
      geom_ribbon(aes(x = time, ymin = m$harvest - err$harvest, 
                      ymax = m$harvest + err$harvest),
                  fill = "darkgreen", alpha = 0.4)  +
      geom_line(aes(time, m$harvest), col = "lightgreen")  #+
#     geom_abline(intercept = e_star, slope = 0, col = "lightgreen", lwd=1,lty=2) # Bifur

## Count how many crashed and add it in a plot title
optimal_crashed = subset(dat, variable == "fishstock" & time == OptTime-1 & value < xT)
p1 <- p1 + opts(title = sprintf("Optimal Harvest dynamics, %d populations crash",
                                dim(optimal_crashed)[1]))
#print(p1)



## Update the main p1 plot to 
## make the crashed trajectories stand out?
#p1 <- p1 + geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "harvest" & 
#          (L1 %in% optimal_crashed$L1)),  col = "darkgreen", alpha = 0.5) +
#     geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "fishstock" &
#          (L1 %in% optimal_crashed$L1)),  col = "darkblue", alpha = 0.5) 



crashed = subset(dat, variable =="unharvested" & time == OptTime - 1 & value < xT)
p2 <- ggplot(data = subset(dat, variable == "unharvested"),
             aes(time, value, group = L1)) + geom_line(alpha = 0.2) + 
  # shows the mean +/- mult * sd , requires Hmisc
  stat_summary(mapping = aes(group = 1), fun.data = mean_sdl, 
               geom = "smooth", mult = 1) +
  opts(title=sprintf("Unfished dynamics, %d populations crash", dim(crashed)[1]))
p2

## Profits plot

p3 <- ggplot(subset(dat, variable == "harvest"), aes(time, profit(value,K), 
                                                     group=L1)) + 
  geom_line(alpha=.2) + labs(x="Time (yrs)", y="Profit" ) + 
  stat_summary(fun.data = mean_sdl, geom="smooth", mapping=aes(group = 1),
               lwd=1, col="darkred", mult=1)
# fun.data should be used for functions that give mean+sd back
# otherwise, use fun.y = mean, fun.ymin = 


cash <- cast(subset(dat,variable=="harvest"), time ~ variable, profit, K) # profits for each rep, by timestep

# add mean line by hand
# p3 <- p3+geom_line(aes(time,rowMeans(cash))) # average profit made as function of time


p4 <- qplot(colSums(cash), xlab="Total Profit", ylab=NULL) + # histogram of total profit made
  geom_vline(xintercept=mean(colSums(cash))) + # expected total profit
  opts(plot.margin=unit(rep(0,4), "lines")) + theme_gray(9) #appearnces
subvp <- viewport(width=.3, height=.3, x=.8, y=.8)

png("profit.png")
print(p3)
print(p4, vp=subvp)
dev.off()
ggsave("samplerun.png", plot=p0)
ggsave("fished.png", plot=p1)
ggsave("unfished.png", plot=p2)

require(socialR)
upload("profit.png samplerun.png fished.png unfished.png", script="Reed.R", tag="PDG_Control")
