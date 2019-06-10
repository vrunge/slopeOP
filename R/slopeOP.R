
#' slopeOP
#' @description Optimal partitioning algorithm for change-in-slope problem with a finite number of states (initial and final values of each segment is restricted to a finite set of values).
#' The algorithm takes into account a continuity constraint between successive segments
#' @param data vector of data to segment
#' @param states vector of states = accessible ending values
#' @param penalty the penalty coefficient. A positive number
#' @param type string defining the pruning type to use
#' @return a list of two vectors and a value  = (changepoints, state parameters, global cost)
#' 'changepoints' is the vector of changepoints (we give the last element of each segment (+ 1 the starting index))
#' 'states' is the vector of successive final values of each segment (+ the very first value)
#' 'globalCost' is a number equal to the global cost of the penalized changepoint problem
slopeOP <- function(data = c(0), states = c(0), penalty = 0, type = "null")
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(data)){stop('data values are not all numeric')}
  if(!is.numeric(states)){stop('states are not all numeric')}
  if(is.unsorted(states)){stop('states should be an increasing vector')}
  if(type != "null" && type != "channel"){stop('Arugment "type" not appropriate. Choose among "null", "channel"')}
  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(penalty < 0){stop('penalty must be nonnegative')}

  ###CALL Rcpp functions###
  res <- slopeOPtransfer(data, states, penalty, type)

  ###Response class slopeOP###
  response <- list(changepoints = res$changepoints + 1, parameters = res$parameters, globalCost = res$globalCost)
  attr(response, "class") <- "slopeOP"

  return(response)
}
