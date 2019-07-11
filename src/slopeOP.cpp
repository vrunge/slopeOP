// MIT License
// Copyright (c) 2019 Vincent Runge

#include<Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include<math.h>
#include"Omega.h"
#include"peltcc_template.h" //Marco Pascucci code
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List slopeOPtransfer(std::vector<double> data, std::vector<double> states, double penalty, std::string constraint = "null", double minAngle = 0, std::string type = "channel")
{
  Omega omega = Omega(states, penalty, data.size());
  //DIFFERENT PRUNING
  if(type == "null" && constraint == "null"){omega.algo(data);}
  if(type == "channel" && constraint == "null"){omega.algoChannel(data);}
  if(type == "pruning" && constraint == "null"){omega.algoPruning(data);}
  if(type == "pruning2" && constraint == "null"){omega.algoPruning2(data);}

  //DIFFERENT CONSTRAINTS
  if(constraint == "up"){omega.algoChannelUP(data);}
  if(constraint == "updown"){omega.algoUPDOWM(data);}
  if(constraint == "smoothing"){omega.algoSMOOTHING(data, minAngle);}

  omega.backtracking(data.size());

  /// RETURN
  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["parameters"] = omega.GetParameters(),
    _["globalCost"] = omega.GetGlobalCost(),
    _["pruningPower"] = omega.GetPruning()
  );
  return res;
}


// MIT License
// Copyright (c) 2019 Marco Pascucci

// [[Rcpp::export]]
List linearOP(std::vector<double> x, std::vector<double> data, double penalty, bool cc = false)
{
  if (x.size() != data.size()) {
    stop("x and y must have the same length.");
  }

  PeltResult<double,double> pr;
  if(cc == false){
    pr = pelt(x,data,penalty);
  } else {
    pr = peltcc(x,data,penalty);
  }

  for(unsigned int i = 0; i < pr.cp.size(); i++){
    pr.cp[i] = pr.cp[i] + 1;
  }

  // RETURN
  List res = List::create(
    _["cp_indexes"] = pr.cp,
    _["x"] = pr.x,
    _["y"] = pr.y,
    _["globalCost"] = pr.cost
  );

  return(res);
}
