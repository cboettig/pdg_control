\name{BevHolt}
\alias{BevHolt}
\title{Harvested Beverton Holt growth model}
\usage{
  BevHolt(x, h, p)
}
\arguments{
  \item{x}{fish population}

  \item{h}{havest level}

  \item{p}{parameters of the growth function, c(A, B),
  where A is the maximum growth rate and B the half-maximum
  in Beverton-Holt.}
}
\value{
  population next year
}
\description{
  Harvested Beverton Holt growth model
}
\details{
  Harvesting takes place before reproduction in ths model.
  The carrying capacity is K <- (pars[1]-1)/pars[2] Try
  with pars <- c(2,4)
}

