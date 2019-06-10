#ifndef COSTS_H
#define COSTS_H

#include <math.h>

class Costs
{
public:
  Costs();
  double slopeCost(double& u, double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& S2t, double& S2T, double& SPt, double& SPT);
  double vhat(double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& SPt, double& SPT);
  unsigned int closestState(double& v, double* states, unsigned int p);
};
#endif // COSTS_H
