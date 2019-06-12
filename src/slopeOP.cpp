#include<Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include<math.h>
#include"Omega.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List slopeOPtransfer(std::vector<double> data, std::vector<double> states, double penalty, std::string constraint = "null", double minAngle = 0, std::string type = "channel")
{
  Omega omega = Omega(states, penalty, data.size());
  if(type == "null" && constraint == "null"){omega.algo(data);}
  if(type == "channel" && constraint == "null"){omega.algoChannel(data);}
  if(type == "pruning" && constraint == "null"){omega.algoPruning(data);}
  if(constraint == "up"){omega.algoChannelUP(data);}
  if(constraint == "updown"){omega.algoUPDOWM(data);}
  if(constraint == "smoothing"){omega.algoSMOOTHING(data, minAngle);}

  omega.backtracking(data.size());

  /// RETURN
  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["parameters"] = omega.GetParameters(),
    _["globalCost"] = omega.GetGlobalCost()
  );
  return res;
}
