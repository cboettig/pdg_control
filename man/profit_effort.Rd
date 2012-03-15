\name{profit_effort}
\alias{profit_effort}
\title{Define a profit function, price minus cost}
\usage{
  profit_effort(price = 1, c0 = 0.1, c1 = 0)
}
\arguments{
  \item{price}{market price (Note, optimal will scrap xT if
  price is high enough!)}

  \item{cost_stock_effect}{fishing extraction costs (per
  unit effort)}

  \item{operating_cost}{costs of having a given harvest
  level (must be less than price)}
}
\value{
  a function computing profit at effort intensity h_i over
  possible stock values x_grid, profit(x_grid, h_i)
}
\description{
  Define a profit function, price minus cost
}

