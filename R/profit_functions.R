#' Define a profit function, price minus cost
#' @param price_fish market price (Note, optimal will scrap xT if price is
#' high enough!) 
#' @param cost_stock_effect fishing extraction costs (per unit effort)
#' @param operating_cost costs of having a given harvest level (must
#' be less than price) 
##' @return a function computing profit at effort intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_effort <-  function(price = 1, c0 = .1, c1 = 0.0){
    function(x_grid, e_i){ 
      sapply(x_grid, function(x_i){
        price * min(e_i * x_i, x_i) -  (c0 + c1 * e_i ) * e_i
    })
  }
}

#' Define a profit function, price minus cost
#' @param price_fish market price (Note, optimal will scrap xT if price is
#' high enough!) 
#' @param c0 Cost, linear with effort (cost = (c0 + c1*E)*E
#' @param c1 Cost, quadratic with effort (cost = (c0 + c1*E)*E
#' @return a function computing profit at harvest intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_harvest  <- function(price = 1, c0 = .1, c1 = 0.0){
  function(x_grid, h_i){
    sapply(x_grid, function(x_i){
      effort <- h_i / (x_i + 1e-12) 
      price_fish * min(h_i, x_i) -  (c0 + c1 * effort ) * effort
      # 1e-12 to avoid NaNs at zero stock condition
    })
  }
}



