// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef OMEGA_H
#define OMEGA_H

#include <math.h>
#include <vector>
#include <list>

#ifndef BUILDING_PYTHON_MODULE
#include "Rcpp.h"
#endif//BUILDING_PYTHON_MODULE

#include "math.h"

class Omega
{
  public:
    Omega(std::vector< double >& values, double beta, unsigned int n);
    ~Omega();

    std::vector< int > GetChangepoints() const;
    std::vector< double > GetParameters() const;
    double GetGlobalCost() const;
    double GetPruning() const;

    void algo(std::vector< double >& data);
    void algoChannel(std::vector< double >& data);
    void algoPruning(std::vector< double >& data);
    void algoPruningMyList(std::vector< double >& data);
    void backtracking(unsigned int n);

    ///////////////
    void algoISOTONIC(std::vector< double >& data);
    void algoUNIMODAL(std::vector< double >& data);
    void algoSMOOTHING(std::vector< double >& data, double minAngle);

  private:
    double penalty;
    unsigned int nbStates;
    double* states;

    double** Q;
    unsigned int** lastIndState;
    unsigned int** lastChpt;

    std::vector< int > changepoints; ///vector of changepoints build by fpop (first index of each segment). size c
    std::vector< double > parameters; ///vector of means build by fpop. size c
    double globalCost;
    double pruning; /// between 0 and 1. 1 if no pruning.
};

#endif // OMEGA_H
