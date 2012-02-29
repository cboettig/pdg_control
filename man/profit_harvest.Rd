\name{profit_harvest}
\alias{profit_harvest}
\title{Define a profit function, price minus cost}
\usage{
  profit_harvest(p = 1, c = 0.001)
}
\arguments{
  \item{p}{price of fish (Note, optimal will scrap xT if
  price is high enough!)}

  \item{c}{fishing extraction costs (per unit effort)}
}
\value{
  a function computing profit at harvest intensity h_i over
  possible stock values x_grid, profit(x_grid, h_i)
}
\description{
  Define a profit function, price minus cost
}

