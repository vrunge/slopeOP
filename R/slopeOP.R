
#' slopeOP
#' @description Optimal partitioning algorithm for change-in-slope problem with a finite number of states (initial and final values of each segment is restricted to a finite set of values).
#' The algorithm takes into account a continuity constraint between successive segments and thus infers a continuous piecewise linear signal
#' @param data vector of data to segment
#' @param states vector of states = set of accessible starting/ending values for segments
#' @param penalty the penalty value A positive number
#' @param constraint string defining a constraint : "null", "up", "updown" or "smoothing"
#' @param minAngle a minimal inner angle in degree between consecutive segments in case constraint = "smoothing"
#' @param type string defining the pruning type to use. "null" = no pruning, "channel" = use monotonicity property or "pruning"
#' @return a list of three elements  = (changepoints, state parameters, global cost)
#' 'changepoints' is the vector of changepoints (we give the extremal values of all segments from left to right)
#' 'states' is the vector of successive states. states[i] is the value we infered at position changepoints[i]
#' 'globalCost' is a number equal to the global cost of the penalized changepoint problem
slopeOP <- function(data = c(0), states = c(0), penalty = 0, constraint = "null", minAngle = 0, type = "channel")
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(data)){stop('data values are not all numeric')}
  if(!is.numeric(states)){stop('states are not all numeric')}
  if(is.unsorted(states)){stop('states should be an increasing vector')}

  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(penalty < 0){stop('penalty must be nonnegative')}

  if(!is.double(minAngle)){stop('minAngle is not a double.')}
  if(minAngle < 0 || minAngle > 180){stop('minAngle must lie between 0 and 180')}

  if(constraint != "null" && constraint != "up" && constraint != "updown" && constraint != "smoothing")
    {stop('Arugment "constraint" not appropriate. Choose among "null", "up", "down" and "smoothing"')}
  if(type != "null" && type != "channel" && type != "pruning" && type != "pruning2")
    {stop('Arugment "type" not appropriate. Choose among "null", "channel" and "pruning"')}



  ###CALL Rcpp functions###
  res <- slopeOPtransfer(data, states, penalty, constraint, minAngle, type)

  ###Response class slopeOP###
  response <- list(changepoints = res$changepoints + 1, parameters = res$parameters, globalCost = res$globalCost)
  attr(response, "class") <- "slopeOP"

  return(response)
}




#' slopeData
#' @description Generate data with a given continuous piecewise linear model
#' @param index a vector of increasing changepoint indexes
#' @param states vector of successive states
#' @param noise noise level = standard deviation of an additional normal noise
#' @return a vector of simulated data
slopeData <- function(index = c(0), states = c(0), noise = 0)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(index)){stop('data values are not all numeric')}
  if(is.unsorted(index)){stop('index should be an increasing vector')}
  if(!is.numeric(states)){stop('states are not all numeric')}
  if(length(index) != length(states)){stop('index and states vectors are of different size')}
  if(!is.double(noise)){stop('noise is not a double.')}
  if(noise < 0){stop('noise must be nonnegative')}

  steps <- diff(states)/diff(index)
  response <- rep(steps,diff(index))
  response <- c(states[1],cumsum(response)+ states[1])
  response <- response + rnorm(length(response),0, noise)

  return(response)
}




