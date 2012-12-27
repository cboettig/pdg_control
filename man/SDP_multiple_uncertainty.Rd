\name{SDP_multiple_uncertainty}
\alias{SDP_multiple_uncertainty}
\title{SDP under multiple uncertainty}
\usage{
  SDP_multiple_uncertainty(f, p, x_grid, h_grid, Tmax = 25,
    sigmas = c(sigma_g = 0.3, sigma_m = 0, sigma_i = 0),
    pdfn = pdfn, profit = function(x, h) pmin(x, h))
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

  \item{sigma}{is the shape parameters for noise
  distribution (sigma_g, sigma_m, sigma_i) (default is no
  noise)}

  \item{pdfn}{is the probability density function (same
  functional form is used for growth, measure, implement).
  (Default is uniform)}

  \item{profit}{is the profit function (defaults to the
  realized harvest)}
}
\value{
  The D matrix giving the optimal harvest for each possible
  state in each timestep. Hence the optimal policy at time
  t is given by the policy function D[,t]. Values of D[,t]
  correspond to the index of h_grid.  Indices of of D[,t]
  correspond to states in y_grid.
}
\description{
  Computes the SDP solution under the case of growth noise,
  implementation errors in harvesting, and meaurement
  errors in the stock assessment.
}

