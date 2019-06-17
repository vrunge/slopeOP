// MIT License
// Copyright (c) 2019 Vincent Runge

#include "Costs.h"
#include<iostream>

#include "math.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Costs::Costs(){}

double Costs::slopeCost(double& u, double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& S2t, double& S2T, double& SPt, double& SPT)
{
  double res = S2T-S2t + (v*v-u*u)/2.0 + (T-t)*(u*u + u*v + v*v)/3.0 + (v-u)*(v-u)/(6.0*(T-t)) - (2.0/(T-t))*(((T+1)*u-(t+1)*v)*(S1T-S1t) + (v-u)*(SPT-SPt));
  return(res);
}


double Costs::vhat(double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& SPt, double& SPT)
{
  double res = (6.0/((T-t-1)*(2.0*(T-t)-1)))*((T+1)*(S1T-S1t) - (SPT-SPt)) - v*(T-t+1)/(2.0*(T-t)-1);
  return(res);
}


unsigned int Costs::closestState(double& v, double* states, unsigned int p)
{
  unsigned int index = p;
  if(v <= states[0]){index = 0;}else{if(v >= states[p-1]){index = p-1;}}
  if(index == p)
  {
    index = 0;
    while(v > states[index]){index = index + 1;}
    if(states[index] + states[index-1] > 2*v){index = index - 1;}
  }
  return(index);
}


bool Costs::pruningTest(unsigned int& tau, unsigned int& t, unsigned int& testT, double& delta, double& DELTA, double& K, double& v)
{
  bool response = false;
  double res = (delta - (DELTA/3.0))*(testT+1) - (v + delta)*(t+2) + (v + (DELTA/3.0))*(tau+1) + ((2*K)+ (DELTA/6.0))/(t-tau);
  if(DELTA > 0 && res <= 0){response = true;}
  if(DELTA < 0 && res >= 0){response = true;}
  if(DELTA == 0){response = false;} /// to be done
  return(response);
}


bool Costs::angleTest(unsigned int& t1, unsigned int& t2, unsigned int& t3, double& v1, double& v2, double& v3, double& minAngle)
{
  bool response = false;

  double angle1 = atan2(v2 - v1, 1.0*(t2 - t1));
  double angle2 = atan2(v3 - v2, 1.0*(t3 - t2));

  double theta = fabs(angle1 - angle2) *  180.0 / M_PI; // in degree
  if(theta <= (180-minAngle)){response = true;}
  //std::cout << " angle1 " << angle1 *  180.0 / M_PI  << " angle2 " << angle2*  180.0 / M_PI << " t1 " << t1 << " v1 " << v1 << " t2 " << t2 << " v2 " << v2 << " t3 " << t3 << " v3 " << v3 << std::endl;

  return(response);
}

