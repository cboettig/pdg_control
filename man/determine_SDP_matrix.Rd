\name{determine_SDP_matrix}
\alias{determine_SDP_matrix}
\title{Determine the Stochastic Dynamic Programming matrix.}
\usage{
  determine_SDP_matrix(f, p, x_grid, h_grid, sigma_g)
}
\arguments{
  \item{f}{the growth function of the escapement population
  (x-h) should be a function of f(t, y, p), with parameters
  p}

  \item{p}{the parameters of the growth function}

  \item{x_grid}{the discrete values allowed for the
  population size, x}

  \item{h_grid}{the discrete values of harvest levels to
  optimize over}

  \item{sigma_g}{the variance of the population growth
  process}
}
\value{
  the transition matrix at each value of h in the grid.
}
\description{
  Determine the Stochastic Dynamic Programming matrix.
}
\details{
  this analytical approach doesn't reliably support other
  sources of variation.  The quality of the analytic
  approximations (lognormal) can be tested.
}

