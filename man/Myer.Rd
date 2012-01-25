\name{Myer}
\alias{Myer}
\title{Discrete-time model with an allee effect for alpha > 1}
\usage{
  Myer(x, h, p)
}
\arguments{
  \item{x}{the current population level}

  \item{h}{harvest effort}

  \item{p}{vector of parameters c(r, alpha, K)}
}
\value{
  the population level in the next timestep
}
\description{
  Discrete-time model with an allee effect for alpha > 1
}
\details{
  A Beverton-Holt style model with Allee effect.  note that
  as written, h is fishing EFFORT, not harvest.  Effort
  above a certain value introduces a fold bifurcation.
  Unharvested carrying capacity is: K <- p[1] * p[3] / 2 +
  sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 The
  (unharvested) allee theshold is given by: x = p[1] * p[3]
  / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
  Bifurcation pt is h = (p[1]*sqrt(p[3])-2)/2 Try with pars
  = c(1,2,6), h=.01

  Consider updating to be a function of x-h, instead?
}

