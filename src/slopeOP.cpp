#include<Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include<math.h>
#include"Omega.h"
using namespace Rcpp;
using namespace std;

//' slopeOP
//' @description Optimal partitioning algorithm for change-in-slope problem with a finite number of states (initial and final values of each segment is restricted to a finite set of values)
//' @param data vector of data to segment
//' @param states vector of states = accessible ending values
//' @param penalty the penalty coefficient. A positive number
//' @param type string defining the pruning type to use
//' @return a list of two vectors and a value  = (changepoints, state parameters, global cost)
//' 'changepoints' is the vector of changepoints (we give the last element of each segment)
//' 'states' is the vector of successive final values of each segment (+ the very first value)
//' 'globalCost' is a number equal to the global cost of the penalized changepoint problem
// [[Rcpp::export]]
List slopeOP(std::vector<double> data, std::vector<double> states, double penalty, std::string type = "null")
{
  ////////////////// //////
  //////// WARNING ////////
  /////// //////////  /////
  if(penalty <= 0){throw std::range_error("Penalty should be a positive number");}
  for(unsigned int i = 0; i < states.size()-1; i++){if(1.0*states[i] >= 1.0*states[i+1]){throw std::range_error("states should be an increasing vector");}}
  if(type != "null" && type != "channel"){throw std::range_error("type is null or channel");}
  //////////
  //////////

  Omega omega = Omega(states, penalty);
  if(type == "null"){omega.algo(data);}
  if(type == "channel"){omega.algoChannel(data);}
  omega.backtracking(data.size());

  /// RETURN
  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["parameters"] = omega.GetParameters(),
    _["globalCost"] = omega.GetGlobalCost()
  );

  return res;
}
