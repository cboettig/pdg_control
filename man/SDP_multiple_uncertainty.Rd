\name{SDP_multiple_uncertainty}
\alias{SDP_multiple_uncertainty}
\title{Determine the SDP matrix for uniform noise under multiple uncertainty}
\usage{
  SDP_multiple_uncertainty(f, p, x_grid, h_grid, sigma_g,
    pdfn = function(P, mu, s) {     if (mu == 0) 
            if (P == 0)              1         else 0
        else dunif(P, mu * (1 - s), mu * (1 + s)) },
    sigma_m, sigma_i, debug = FALSE)
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

  \item{sigma_g}{is the shape parameter (width) of the
  multiplicitive growth noise}

  \item{pdfn}{is the shape of the growth noise, which need
  not be uniform (is by default)}

  \item{sigma_m}{is the half-width of the uniform
  measurement error (assumes uniform distribution)}

  \item{sigma_i}{is the half-width of the implementation
  noise (always assumes uniform distribution)}
}
\description{
  Computes the transition matrix under the case of growth
  noise, implementation errors in harvesting, and
  meaurement errors in the stock assessment.
}

