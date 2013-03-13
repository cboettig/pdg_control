\name{coral}
\alias{coral}
\title{Coral-Parrotfish model}
\usage{
  coral(x, h,
    p = c(a = 0.1, g = 1, T = , gamma = 0.8, r = 1, d = 0.44, s = 0.49, K = 1))
}
\arguments{
  \item{x}{vector of population levels: macroalgae, coral,
  parrotfish}

  \item{h}{harvesting effort on parrotfish}

  \item{p}{c(a, g, T, gamma, r, d, s, K) 1 2 3 4 5 6 7 8}
}
\description{
  Coral-Parrotfish model
}
\details{
  $$\begin{align} \frac{dM}{dt} = aMC - \frac{g(P) M}{M+T}
  + \gamma M T \\ \frac{dC}{dt} = rTC - dC - a M C \\
  \frac{dP}{dt} = sP \left( 1- \frac{P}{\beta K(C) }
  \right) - h P \end{align}$$
}
\references{
  Blackwood et al. (2011) doi:10.1890/10-2195.1
}

