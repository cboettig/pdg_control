\name{BevHolt_effort}
\alias{BevHolt_effort}
\title{Beverton Holt growth model with fishing-effort based control}
\usage{
  BevHolt_effort(x, h, p)
}
\arguments{
  \item{x}{fish population}

  \item{h}{fishing effort (harvest is effort times fish
  population)}

  \item{p}{parameters of the growth function, c(A, B),
  where A is the maximum growth rate and B the half-maximum
  in Beverton-Holt.}
}
\value{
  population next year
}
\description{
  Beverton Holt growth model with fishing-effort based
  control
}
\details{
  Harvesting takes place before reproduction in ths model.
  The carrying capacity is K <- (pars[1]-1)/pars[2] Try
  with pars <- c(2,4)
}

