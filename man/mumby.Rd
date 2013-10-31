\name{mumby}
\alias{mumby}
\title{Coral-Parrotfish model}
\usage{
  mumby(x, h,
    p = c(a = 0.1, g = 0.5, gamma = 0.8, r = 1, d = 0.44),
    dt = 0.025)
}
\arguments{
  \item{x}{vector of population levels: macroalgae, coral,
  parrotfish}

  \item{h}{harvesting effort on parrotfish}

  \item{p}{c(a, g, gamma, r, d)}

  \item{dt}{descrization scale}
}
\description{
  Coral-Parrotfish model
}
\details{
  \deqn{ \frac{dM}{dt} = aMC - \frac{(g-h) M}{M+T} + \gamma
  M T \\ \frac{dC}{dt} = rTC - dC - a M C \\ T = 1 - M - C
  }{}
}
\references{
  Mumby et al. (2007) doi:10.1038/nature06252 Blackwood et
  al. (2011) doi:10.1890/10-2195.1
}

