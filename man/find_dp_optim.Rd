\name{find_dp_optim}
\alias{find_dp_optim}
\title{Identify the dynamic optimum using backward iteration (dynamic programming)}
\usage{
  find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT,
    profit, delta, reward = 10, interval = 1)
}
\arguments{
  \item{SDP_Mat}{the stochastic transition matrix at each h
  value}

  \item{x_grid}{the discrete values allowed for the
  population size, x}

  \item{h_grid}{the discrete values of harvest levels to
  optimize over}

  \item{OptTime}{the stopping time}

  \item{xT}{the boundary condition population size at
  OptTime}

  \item{c}{the cost/profit function, a function of
  harvested level}

  \item{delta}{the exponential discounting rate}

  \item{reward}{the profit for finishing with >= Xt fish at
  the end (i.e. enforces the boundary condition)}

  \item{interval}{the delay interval between updated
  assessment/quotas}
}
\value{
  list containing the matrices D and V.  D is an x_grid by
  OptTime matrix with the indices of h_grid giving the
  optimal h at each value x as the columns, with a column
  for each time.  V is a matrix of x_grid by x_grid, which
  is used to store the value function at each point along
  the grid at each point in time.  The returned V gives the
  value matrix at the first (last) time.
}
\description{
  Identify the dynamic optimum using backward iteration
  (dynamic programming)
}

