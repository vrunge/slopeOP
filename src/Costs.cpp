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
  ///REMARK : t -> t+1 and T -> T+1 to get indexation starting from 1
  double res = S2T-S2t + (v*v-u*u)/2.0 + (T-t)*(u*u + u*v + v*v)/3.0 + (v-u)*(v-u)/(6.0*(T-t)) - (2.0/(T-t))*(((T+1)*u-(t+1)*v)*(S1T-S1t) + (v-u)*(SPT-SPt));
  return(res);
}


//####### vhat #######////####### vhat #######////####### vhat #######//
//####### vhat #######////####### vhat #######////####### vhat #######//

double Costs::vhat(double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& SPt, double& SPT)
{
  ///REMARK : t -> t+1 and T -> T+1 to get indexation starting from 1
  double res = (6.0/((T-t-1)*(2.0*(T-t)-1)))*((T+1)*(S1T-S1t) - (SPT-SPt)) - v*(T-t+1)/(2.0*(T-t)-1);
  return(res);
}


//####### closestStateIndex #######////####### closestStateIndex #######////####### closestStateIndex #######//
//####### closestStateIndex #######////####### closestStateIndex #######////####### closestStateIndex #######//

unsigned int Costs::closestStateIndex(double& v, double* states, unsigned int p)
{
  if(v <= states[0]){return(0);}
  if(v >= states[p - 1]){return(p - 1);}

  // binary search
  unsigned int i = 0, j = p, mid = 0;
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

bool Costs::pruningTest(unsigned int& tau, unsigned int& t, unsigned int& testT, double& delta, double& DELTA, double& K, double& v)
{
  bool response = false;
  double res = (delta - (DELTA/3.0))*(testT+1) - (v + delta)*(t+2) + (v + (DELTA/3.0))*(tau+1) + ((2*K)+ (DELTA/6.0))/(t-tau);
  if(DELTA > 0 && res <= 0){response = true;}
  if(DELTA < 0 && res >= 0){response = true;}
  if(DELTA == 0){response = false;} /// to be done
  return(response);
}


//####### angleTest #######////####### angleTest #######////####### angleTest #######//
//####### angleTest #######////####### angleTest #######////####### angleTest #######//

bool Costs::angleTest(unsigned int& t1, unsigned int& t2, unsigned int& t3, double& v1, double& v2, double& v3, double& minAngle)
{
  bool response = false;

  double angle1 = atan2(v2 - v1, 1.0*(t2 - t1));
  double angle2 = atan2(v3 - v2, 1.0*(t3 - t2));

  double theta = fabs(angle1 - angle2) *  180.0 / M_PI; // in degree
  if(theta <= (180 - minAngle)){response = true;}

  return(response);
}
