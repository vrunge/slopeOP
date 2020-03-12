// MIT License
// Copyright (c) 2019 Vincent Runge

#include<iostream>
#include "Costs.h"
#include "math.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Costs::Costs(){}

//####### slopeCost #######////####### slopeCost #######////####### slopeCost #######//
//####### slopeCost #######////####### slopeCost #######////####### slopeCost #######//

double Costs::slopeCost(double& u, double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& S2t, double& S2T, double& SPt, double& SPT)
{
  return(S2T-S2t - (2.0/(T-t))*((T*u-t*v)*(S1T-S1t) + (v-u)*(SPT-SPt)) + (v*v-u*u)/2.0 + (T-t)*(u*u + u*v + v*v)/3.0 + (v-u)*(v-u)/(6.0*(T-t)));
}

//####### vhat #######////####### vhat #######////####### vhat #######//
//####### vhat #######////####### vhat #######////####### vhat #######//

double Costs::vhat(double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& SPt, double& SPT)
{
  return((6.0/((T-t-1)*(2.0*(T-t)-1)))*(T*(S1T-S1t) - (SPT-SPt)) - v*(T-t+1)/(2.0*(T-t)-1));
}

//####### closestStateIndex #######////####### closestStateIndex #######////####### closestStateIndex #######//
//####### closestStateIndex #######////####### closestStateIndex #######////####### closestStateIndex #######//

unsigned int Costs::closestStateIndex(double& v, double* states, unsigned int p)
{
  if(v <= states[0]){return(0);}
  if(v >= states[p - 1]){return(p - 1);}

  // binary search
  unsigned int i = 0;
  unsigned int j = p;
  unsigned int mid = 0;
  while(i < j)
  {
    mid = (i + j)/2;
    if(states[mid] == v){return(mid);}

    if(v < states[mid])
    {
      if(mid > 0 && v > states[mid - 1])
        {if(states[mid - 1] + states[mid] > 2*v){return(mid - 1);}else{return(mid);}}
      j = mid;
    }
    else
    {
      if (mid < (p - 1) && v < states[mid + 1])
        {if(states[mid] + states[mid + 1] > 2*v){return(mid);}else{return(mid + 1);}}
      i = mid + 1;
    }
  }
  return(mid);
}


//####### pruningTest #######////####### pruningTest #######////####### pruningTest #######//
//####### pruningTest #######////####### pruningTest #######////####### pruningTest #######//

bool Costs::pruningTest(double r, double s, unsigned int tau, unsigned int t, unsigned int n, double S, double Ap, double Am, double Gp, double Gm)
{
  bool response = false;
  if(r == s){response = true;}
  else if(s > r)
  {
    double state1 = (r+2*s)/6.0;
    double tempGp = tau*state1 - (s-r)/(12.0*(t-tau)) + S/(1.0*(t-tau)) + Gp;
    double ApC = Ap - state1;
    if((ApC*(t+1) + tempGp > 0) && (ApC*n + tempGp > 0)){response = true;}
  }
  else /// s < r
  {
    double state2 = (r+2*s)/6.0;
    double tempGm = tau*state2 - (s-r)/(12.0*(t-tau)) + S/(1.0*(t-tau)) + Gm;
    double AmC = Am - state2;
    if((AmC*(t+1) + tempGm > 0) && (AmC*n + tempGm > 0)){response = true;}
  }
  return(response);
}


//####### angleTest #######////####### angleTest #######////####### angleTest #######//
//####### angleTest #######////####### angleTest #######////####### angleTest #######//

bool Costs::angleTest(unsigned int& t1, unsigned int& t2, unsigned int& t3, double& v1, double& v2, double& v3, double& minAngle)
{
  bool response = false;
  double cosAngleRad = ((1.0*t1-1.0*t2)*(1.0*t3-1.0*t2) + (v1-v2)*(v3-v2))/(sqrt(((1.0*t1-1.0*t2)*(1.0*t1-1.0*t2) + (v1-v2)*(v1-v2))*((1.0*t3-1.0*t2)*(1.0*t3-1.0*t2) + (v3-v2)*(v3-v2))));
  double theta = acos(cosAngleRad) *180.0 / M_PI; // in degree

  if(theta >= minAngle){response = true;}
  if((t1 == t2) && (v1 == v2)){response = true;}
  //std::cout << "t1 " << t1 << " v1 "<< v1 << " t2 "<< t2 << " v2 "<< v2 << " t3 " << t3 << " v3 "<< v3 << " response " << response << " && ";

  return(response);
}
