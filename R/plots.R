#' Show the ensemble fishstock and harvest dynamics
#' @param sims simulation results
#' @param S Reed's S, from opt$S, where opt <- find_dp_optim(...)
#' @param allee the allee threshold (defaults to 0)
#' @param fish_reps logical - show replicate trajectories of fish stocks?
#' @param fish_sd logical - show ribbon of SD of fishstock?
#' @param fish_mean logical - show mean fishstock?
#' @param harvest - include harvest level on the same plot?
#' @param harvest_reps - show individual reps for harvest?
#' @param harvest_mean - show mean harvest level over time?
#' @param harvest_sd - show harvest SD
#' @param crash_count - count the number of reps that have fallen below allee
#' @return a ggplot object (display with print(object)) 
#' @details if line values for any S, allee, e_star are not specified 
#' then these are omitted from the plot.
#' @import ggplot2
#' @import reshape2
#' @import Hmisc
#' @export 
plot_replicates <- function(sims, S = NULL, allee = 0, e_star = NULL,
                            fish_reps = TRUE, fish_sd = TRUE, 
                            fish_mean = TRUE, harvest = TRUE, 
                            harvest_reps = TRUE, harvest_mean = TRUE,
                            harvest_sd = TRUE, crash_count=TRUE){

  require(ggplot2)

  ## Should rewrite this to be act like a proper ggplot2 extension 

  dat <- melt(sims, id="time") # reshapes the data matrix to "long" form
  dat <<- dat
  # some stats on the replicates, (stat_sumary can do this instead)
  m <- dcast(dat, time ~ variable, mean) # mean population
  m <<- m
  err <- dcast(dat, time ~ variable, sd) # sd population
  err <<- err
  
  ## Create the ggplot object with this data
  p1 <- ggplot(dat) 

  ## Add the layers requested 
  if(fish_reps) # Replicate harvested dynamics
    p1 <- p1+ geom_line(aes(time, value, group = L1), data = 
                  subset(dat, variable == "fishstock"), alpha = 0.2)
  if(fish_sd) 
    p1 <- p1 + geom_ribbon(aes(x = time, 
                           ymin = m$fishstock - err$fishstock,
                           ymax = m$fishstock + err$fishstock),
                           fill = "darkblue", alpha = 0.4)
  if(fish_mean)
    p1 <- p1 + geom_line(aes(time, m$fishstock), col = "lightblue") 
  ## And the same kind of plots for harvest-effort levels
  if(harvest){
    if(harvest_reps)
      p1 <- p1 +geom_line(aes(time, value, group = L1), 
              data = subset(dat, variable == "harvest"), 
                            alpha=.2, col = "darkgreen") 
    if(harvest_sd)
      p1 <- p1 + geom_ribbon(aes(x = time, 
                             ymin = m$harvest - err$harvest, 
                             ymax = m$harvest + err$harvest),
                             fill = "darkgreen", alpha = 0.4)
    if(harvest_mean)
      p1 <- p1 + geom_line(aes(time, m$harvest), col = "lightgreen")
  }
  if(crash_count){  # Count how many crashed and add it in a plot title
    optimal_crashed = subset(dat, variable == "fishstock" 
                             & time == max(dat["time"]) - 1 & value <= allee)
    p1 <- p1 + opts(title = 
          sprintf("Optimal Harvest dynamics, %d populations crash",
                  dim(optimal_crashed)[1]))
  }
  if(!is.null(S)) # show Reed's S: optimal escapement 
    p1 <- p1 + geom_abline(intercept=S, slope = 0)          
  if(!is.null(allee)) # show Allee threshold
    p1 <- p1 + geom_abline(intercept=allee, slope = 0, lty=2) 
  if(!is.null(e_star)) # Bifurcation level of harvesting 
    geom_abline(intercept = e_star, slope = 0, col = "lightgreen", lwd=1,lty=2)
  # return the plotted object
  p1
}


