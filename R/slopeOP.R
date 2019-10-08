
#' slopeOP
#' @description Optimal partitioning algorithm for change-in-slope problem with a finite number of states (beginning and ending values of each segment is restricted to a finite set of values).
#' The algorithm takes into account a continuity constraint between successive segments and infers a continuous piecewise linear signal.
#' @param data vector of data to segment
#' @param states vector of states = set of accessible starting/ending values for segments in increasing order.
#' @param penalty the penalty value (a positive number)
#' @param constraint string defining a constraint : "null", "isotonic", "unimodal" or "smoothing"
#' @param minAngle a minimal inner angle in degree between consecutive segments in case constraint = "smoothing"
#' @param type string defining the pruning type to use. "null" = no pruning, "channel" = use monotonicity property or "pruning" = pelt-type property
#' @param testMode a boolean, if true the function also returns the percent of elements to scan (= ratio scanned elements vs. scanned elements if no pruning)
#' @return a list of three elements  = (changepoints, state parameters, global cost)
#' \itemize{
#'   \item \strong{changepoints} is the vector of changepoints (we give the extremal values of all segments from left to right)
#'   \item \strong{states} is the vector of successive states. states[i] is the value we infered at position changepoints[i]
#'   \item \strong{globalCost} is a number equal to the global cost of the penalized change-in-slope problem
#'   \item \strong{pruning} is the percent of positions to consider in matrix Q  (returned only if testMode = TRUE)
#' }
slopeOP <- function(data = c(0), states = c(0), penalty = 0, constraint = "null", minAngle = 0, type = "channel", testMode = FALSE)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(data)){stop('data values are not all numeric')}
  if(!is.numeric(states)){stop('states are not all numeric')}
  if(is.unsorted(states)){stop('states should in increasing order')}
  if(length(unique(states)) < length(states)){stop('states is not a strictly increasing sequence')}

  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(penalty < 0){stop('penalty must be nonnegative')}
  if(!is.double(minAngle)){stop('minAngle is not a double.')}
  if(minAngle < 0 || minAngle > 180){stop('minAngle must lie between 0 and 180')}

  if(constraint != "null" && constraint != "isotonic" && constraint != "unimodal" && constraint != "smoothing")
    {stop('Arugment "constraint" not appropriate. Choose among "null", "isotonic", "unimodal" and "smoothing"')}
  if(type != "null" && type != "channel" && type != "pruning" && type != "pruningMyList")
    {stop('Arugment "type" not appropriate. Choose among "null", "channel" and "pruning"')}

  if(!is.logical(testMode)){stop('testMode must be a boolean')}

  ###CALL Rcpp functions###
  res <- slopeOPtransfer(data, states, penalty, constraint, minAngle, type)

  ###Response class slopeOP###
  ### ATTENTION : we here remove one penalty to globalCost
  if(testMode == FALSE){response <- list(changepoints = res$changepoints + 1, parameters = res$parameters, globalCost = res$globalCost - penalty)}
  if(testMode == TRUE){response <- list(changepoints = res$changepoints + 1, parameters = res$parameters, globalCost = res$globalCost - penalty, pruning = res$pruningPower)}

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


#' plot.slopeOP
#' @description Plot the result of the slopeOP function with the data
#' @param x a slopeOP class object
#' @param ... Other parameters
#' @param data the data from which we get res
#' @param chpt vector of changepoints of the model
#' @param states vector of states of the model
#' @return plot data and the inferred slopeOP result
plot.slopeOP <- function(x, ..., data, chpt = NULL, states = NULL)
{
  n <- 1:length(data)
  p <- length(x$changepoints)
  xbis <- x$changepoints
  y <- x$parameters

  plot(1:length(data), data, pch = '+')
  for(i in 1:(p-1))
  {
    segments(xbis[i], y[i], xbis[i+1], y[i+1], col= 2, lty = 1, lwd = 5)
  }
  if(length(chpt) > 0 && length(chpt) == length(states))
  {
    q <- length(chpt)
    for(i in 1:(q-1))
    {
      segments(chpt[i], states[i], chpt[i+1], states[i+1], col= 4, lty = 1, lwd = 5)
    }
  }
}

