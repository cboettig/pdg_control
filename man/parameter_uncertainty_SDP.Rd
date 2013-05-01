\name{parameter_uncertainty_SDP}
\alias{parameter_uncertainty_SDP}
\title{Computes the transition matrix under the case of parametric uncertainty}
\usage{
  parameter_uncertainty_SDP(f, p, x_grid, h_grid, pardist,
    n_mc = 100)
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

  \item{pardist}{a matrix with columns as variables and
  rows as the Monte Carlo samples from the posterior}
}
\description{
  Computes the transition matrix under the case of
  parametric uncertainty
}

