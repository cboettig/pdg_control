<!--roptions dev="png", fig.width=7, fig.height=5, fig.path='ex-out-', tidy=FALSE, warning=FALSE, comment=NA, message=FALSE, cache=FALSE-->

<!--begin.rcode echo=FALSE 
render_gfm()
opts_knit$set(upload = TRUE)
## use flickr to upload with these options
require(socialR)
options(flickrOptions=list(
  description="https://github.com/cboettig/pdg_control/blob/master/inst/examples/",
  tags="stochpop, pdg_control"))
opts_knit$set(upload.fun = flickr.url)
end.rcode-->


<!--begin.rcode pars
require(pdgControl)
p_grid = seq(0.01,.99, length=5) 
x_grid = seq(1,10,length=10) 
sigma_g = 0.2
h_grid <- x_grid-1 
T <- 5
xT <- 0
z_g <- function() rlnorm(1,  0, sigma_g) 
profit <- profit_harvest(price=10, c0=30) 
delta <- 0.05
reward <- 0
end.rcode-->

Static solution
<!--begin.rcode static
bevholt <- function(x,h, p) p[1] * (x-h) / (1 - p[2] * (x-h))
sdp <- determine_SDP_matrix(bevholt, c(1.5, 0.05), x_grid, h_grid, .2)
static <- find_dp_optim(sdp, x_grid, h_grid, T, xT=0, profit, delta, reward)
end.rcode-->

Active Adaptive Mangement solution
<!--begin.rcode active
f1 <- setmodel(bevholt, c(1.5,0.05))
f2 <- setmodel(bevholt, c(3,0.05))
M <- model_uncertainty(f1, f2, x_grid, p_grid, h_grid)
active <- dp_optim(M, x_grid, h_grid, T, xT=0, profit, delta, reward, p_grid=p_grid) 

sim <- active_adaptive_simulate(bevholt, c(1.5, 0.05), x_grid, h_grid, p_grid, 
                                x_grid[6], p_grid[4], active$D, z_g, update_belief(f1,f2))
end.rcode-->



