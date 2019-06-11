#include<Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include<math.h>
#include"Omega.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List slopeOPtransfer(std::vector<double> data, std::vector<double> states, double penalty, std::string type = "null")
{
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
