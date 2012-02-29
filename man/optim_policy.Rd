\name{optim_policy}
\alias{optim_policy}
\title{Optimum policy when changing harvest level has an adjustment cost}
\usage{
  optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT,
    profit, delta, reward = 0, P = 0,
    penalty = c("L1", "L2", "asymmetric", "fixed", "none"))
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

  \item{P}{the cost of adjusting a policy, proportional to
  the amount of change}

  \item{penalty}{the kind of penalty applied: currently L1,
  L2, or assymetric}
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

