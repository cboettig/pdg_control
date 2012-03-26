\name{active_adaptive_simulate}
\alias{active_adaptive_simulate}
\title{Forward simulate given the optimal havesting policy, D}
\usage{
  active_adaptive_simulate(f, pars, x_grid, h_grid, p_grid,
    x0, p0, D, z_g, update_belief)
}
\arguments{
  \item{f}{the true growth function of the escapement
  population (x-h) should be a function of f(y, p), with
  parameters p}

  \item{pars}{the parameters of the growth function}

  \item{x_grid}{the discrete values allowed for the
  population size, x}

  \item{h_grid}{the discrete values of harvest levels to
  optimize over}

  \item{p_grid}{the discrete values of belief allowed}

  \item{x0}{initial stock size}

  \item{p0}{initial belief}

  \item{D}{the optimal solution indices on h_grid, given
  for each possible state at each timestep}

  \item{z_g}{a function which returns a random multiple for
  population growth for the implementation uncertainty in
  quotas}

  \item{update_belief}{a function of x[t], p[t], x[t+1]}
}
\value{
  a data frame with the time, fishstock, harvested amount,
  and what the escapement ("unharvested").
}
\description{
  Forward simulate given the optimal havesting policy, D
}

