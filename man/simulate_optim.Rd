\name{simulate_optim}
\alias{simulate_optim}
\title{Forward simulate given the optimal havesting policy when cost depends on prev harvest}
\usage{
  simulate_optim(f, pars, x_grid, h_grid, x0, D, z_g, z_m,
    z_i, alt_D = NULL, OptTime = dim(D[[1]])[2])
}
\arguments{
  \item{f}{the growth function of the escapement population
  (x-h) should be a function of f(y, p), with parameters p}

  \item{p}{the parameters of the growth function}

  \item{x_grid}{the discrete values allowed for the
  population size, x}

  \item{h_grid}{the discrete values of harvest levels to
  optimize over}

  \item{Xo}{initial stock size}

  \item{D}{the optimal solution indices on h_grid, given
  for each possible state at each timestep}

  \item{z_g}{a function which returns a random multiple for
  population growth}

  \item{z_m}{a function which returns the measurement
  uncertainty in the assessment of stock size}

  \item{z_i}{a function which returns a random number from
  a distribution for the implementation uncertainty in
  quotas}

  \item{alt_D}{an option to compute harvest with an
  alternate "optimal" solution as well, for comparison
  (indep of h_prev)}
}
\value{
  a data frame with the time, fishstock, harvested amount,
  and what the escapement ("unharvested") is. Also the
  alternate fishstock and harvest dynamics, if calculated
}
\description{
  This simulates a process where there is a cost to
  changing the policy.
}

