#' Define a profit function, price minus cost
#' @param price_fish market price (Note, optimal will scrap xT if price is
#' high enough!) 
#' @param cost_stock_effect fishing extraction costs (per unit effort)
#' @param operating_cost costs of having a given harvest level (must
#' be less than price) 
##' @return a function computing profit at effort intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_effort <-  function(price_fish = 1, cost_stock_effect = 0.000,
                            operating_cost = 0.1 * price_fish){
#' @param x_grid is a the grid of state values (profit will evaluate at each of them)
#' @param h_i is the current harvesting *effort* (effort*stocksize = catch) 
#' @return the profits of fishing at intensity h_i given the stock value equals x_i
#' for each x_i in the grid.   
#' @details Due to the symmetry, you can actually compute the profit over a range 
#'  of harvest values, rather than a range of stock values, by simply swapping 
#'  x and h, i.e. give a vector of h values as x_grid, and a single stock size as h_i. 
    function(x_grid, h_i){ 
      sapply(x_grid, function(x_i){
        h <- h_i * x_i # harvest propotional to effort 
        price_fish * min(h, x_i) - operating_cost * h  - cost_stock_effect / (x_i + 1e-12)
    })
  }
}

#' Define a profit function, price minus cost
#' @param price_fish market price (Note, optimal will scrap xT if price is
#' high enough!) 
#' @param cost_stock_effect fishing extraction costs (per unit effort)
#' @param operating_cost costs of having a given harvest level (must
#' be less than price) 
#' @return a function computing profit at harvest intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_harvest  <- function(price_fish = 1, cost_stock_effect = 0.000,
                            operating_cost = 0.1 * price_fish){
#' @param x_grid is a the grid of state values (profit will evaluate at each of them)
#' @param h_i is total harvest level
#' @return the profits of harvesting at intensity h_i for each possible stock
#' size in x_grid.  
#' @details Due to the symmetry, you can actually compute the profit over a range 
#'  of harvest values, rather than a range of stock values, by simply swapping 
#'  x and h, i.e. give a vector of h values as x_grid, and a single stock size as h_i. 
  function(x_grid, h_i){
    sapply(x_grid, function(x_i){
      price_fish * min(h_i, x_i) - operating_cost * h_i  -
        cost_stock_effect / (x_i + 1e-12)
      # 1e-12 to avoid NaNs at zero stock condition
    })
  }
}



