<!--roptions dev='png', fig.width=7, fig.height=5, tidy=FALSE, warning=FALSE, message=FALSE, comment=NA, cache.path="policycost/", cache=FALSE-->

<!--begin.rcode setup, include=FALSE
render_gfm()  
opts_knit$set(upload = TRUE)   
require(socialR)
options(flickrOptions=list(
  description="https://github.com/cboettig/pdg_control/blob/master/inst/examples/profile_optimpolicy.md",
  tags="stochpop, pdg_control"))
opts_knit$set(upload.fun = flickr.url)
read_chunk("Reed.R")
end.rcode-->

<!--begin.rcode libraries, echo=FALSE
end.rcode-->



<!--begin.rcode parameters
delta <- 0.01     # SMALLER economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 0       # bonus for satisfying the boundary condition
end.rcode-->
<!--begin.rcode noise_dists
end.rcode-->
<!--begin.rcode BevHolt
end.rcode-->
<!--begin.rcode initx
end.rcode-->
and we use a harvest-based profit function with default parameters
<!--begin.rcode profit_
profit <- profit_harvest(price_fish = 1, cost_stock_effect = 0,
 operating_cost = 0.1)
end.rcode-->

Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 
<!--begin.rcode create_grid
end.rcode-->
<!--begin.rcode determine_SDP_matrix
end.rcode-->

<!--begin.rcode
require(profr)
end.rcode-->

<!--begin.rcode find_dp_optim_ 
fast <- profr(opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward))
end.rcode-->

<!--begin.rcode policycost_optim
slow <- profr(policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, P = .3, penalty = "asym"))
end.rcode-->


<!--begin.rcode
plot(fast)
plot(slow)
fast
slow
end.rcode-->

