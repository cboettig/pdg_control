#' Define a profit function, price minus cost
#' @param p price of fish (Note, optimal will scrap xT if price is high enough!) 
#' @param c fishing extraction costs (per unit effort)
#' @return a function computing profit at effort intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_effort <- function(p = 1, c = 0.001){
#' @param x_grid is a the grid of state values (profit will evaluate at each of them)
#' @param h_i is the current harvesting *effort* (effort*stocksize = catch) 
#' @return the profits of fishing at intensity h_i given the stock value equals x_i
#' for each x_i in the grid.   
#' @details Due to the symmetry, you can actually compute the profit over a range 
#'  of harvest values, rather than a range of stock values, by simply swapping 
#'  x and h, i.e. give a vector of h values as x_grid, and a single stock size as h_i. 
  function(x_grid, h_i){ 
    harvest <- x_grid * h_i 
    sapply(harvest, function(x) max(0, p * x - c / x))
  }
}


#' Define a profit function, price minus cost
#' @param p price of fish (Note, optimal will scrap xT if price is high enough!) 
#' @param c fishing extraction costs (per unit effort)
#' @return a function computing profit at harvest intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_harvest  <- function(p = 1, c = 0.000){
#' @param x_grid is a the grid of state values (profit will evaluate at each of them)
#' @param h_i is total harvest level
#' @return the profits of harvesting at intensity h_i for each possible stock
#' size in x_grid.  
#' @details Due to the symmetry, you can actually compute the profit over a range 
#'  of harvest values, rather than a range of stock values, by simply swapping 
#'  x and h, i.e. give a vector of h values as x_grid, and a single stock size as h_i. 
  function(x_grid, h_i){
    sapply(x_grid, function(x_i){
      max(0, p * min(h_i, x_i) - c / x_i)
    })
  }
}



