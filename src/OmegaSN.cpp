// MIT License
// Copyright (c) 2019 Vincent Runge

#include "OmegaSN.h"
#include "Costs.h"

#include <algorithm> // std::reverse // std::min // std::max
#include<iostream>
#include <stdlib.h>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

OmegaSN::OmegaSN(std::vector< double >& values, unsigned int nbSeg, unsigned int n)
{
  nbStates = values.size();
  states = new double[nbStates];
  for(unsigned int i = 0; i < nbStates; i++){states[i] = values[i];}

  ///
  /// MATRIX INITIALIZATION
  ///
  S12P = new double*[3]; ///matrix of vectors S1, S2 and SP
  Q = new double*[nbStates]; ///matrix of costs
  lastChpt = new unsigned int*[nbStates]; ///matrix of best last changepoints
  lastIndState = new unsigned int*[nbStates]; ///matrix of starting states for the best last segment

  for(unsigned int i = 0; i < 3; i++){S12P[i] = new double[n];}
  for(unsigned int i = 0; i < nbStates; i++){Q[i] = new double[n];}
  for(unsigned int i = 0; i < nbStates; i++){lastChpt[i] = new unsigned int[n];}
  for(unsigned int i = 0; i < nbStates; i++){lastIndState[i] = new unsigned int[n];}

  nbSegments = nbSeg;
}



OmegaSN::~OmegaSN()
{
  delete(states);
  states = NULL;
  for(unsigned int i = 0; i < 3; i++){delete(S12P[i]);}
  for(unsigned int i = 0; i < nbStates; i++){delete(Q[i]);}
  for(unsigned int i = 0; i < nbStates; i++){delete(lastChpt[i]);}
  for(unsigned int i = 0; i < nbStates; i++){delete(lastIndState[i]);}
  delete [] S12P;
  S12P = NULL;
  delete [] Q;
  Q = NULL;
  delete [] lastChpt;
  lastChpt = NULL;
  delete [] lastIndState;
  lastIndState = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > OmegaSN::GetChangepoints() const {return(changepoints);}
std::vector< double > OmegaSN::GetParameters() const {return(parameters);}
double OmegaSN::GetGlobalCost() const {return(globalCost);}
double OmegaSN::GetPruning() const {return(pruning);}


//####### preprocessing #######////####### preprocessing #######////####### preprocessing #######//
//####### preprocessing #######////####### preprocessing #######////####### preprocessing #######//

double** OmegaSN::preprocessing(std::vector< double >& data) const
{
  unsigned int n = data.size();
  S12P[0][0] = data[0];
  S12P[1][0] = data[0] * data[0];
  S12P[2][0] = data[0];
  for(unsigned int i = 1; i < n; i++){S12P[0][i] = S12P[0][i-1] + data[i];}
  for(unsigned int i = 1; i < n; i++){S12P[1][i] = S12P[1][i-1] + (data[i] * data[i]);}
  for(unsigned int i = 1; i < n; i++){S12P[2][i] = S12P[2][i-1] + (i+1) * data[i];}
  return(S12P);
}


// algo // algoChannel // algoPruning // algoPruningPELT // backtracking // algoISOTONIC // algoUNIMODAL // algoOUTLIER

//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//

void OmegaSN::algo(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;
  unsigned int zero = 0;

  /// PREPROCESSING
  S12P = preprocessing(data);

  Costs cost;
  ///
  /// FILL FIRST COLUMN in Q
  ///
  for(unsigned int i = 0; i < p; i++)
  {
    Q[i][0] = (data[0] - states[i])*(data[0] - states[i]);
    lastChpt[i][0] = 0;
    lastIndState[i][0] = 0;
  }

  ///
  /// ALGO
  ///
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 1; T < n; T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*T
      /////
      // INITIALIZATION of temp_Q
      temp_Q = Q[0][0] + cost.slopeCost(states[0], states[v], zero, T, S12P[0][0], S12P[0][T], S12P[1][0], S12P[1][T], S12P[2][0], S12P[2][T]);
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 0; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore column of states
        {
          temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
          if(temp_Q > temp_cost)
          {
            temp_Q = temp_cost;
            temp_indState = u;
            temp_chpt = t;
          }
        }
      }
      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;
    }
  }
  pruning = 1; ///We went through all the elements in matrix Q
}

//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//

void OmegaSN::algoChannel(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;
  unsigned int zero = 0;

  /// PREPROCESSING
  S12P = preprocessing(data);

  Costs cost;
  ///
  /// FILL FIRST COLUMN in Q
  ///
  for(unsigned int i = 0; i < p; i++)
  {
    Q[i][0] = (data[0] - states[i])*(data[0] - states[i]);
    lastIndState[i][0] = 0;
    lastChpt[i][0] = 0;
  }

  ///
  /// ALGO
  ///
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///
  /// CHANNEL INFORMATION
  /// u1 / u2 = "min / max" in each column of Q
  ///
  unsigned int* u1 = new unsigned int[n];
  unsigned int* u2 = new unsigned int[n];
  unsigned int theStart;
  unsigned int theEnd;
  double theV = 0;
  unsigned int indexTheV = 0;

  unsigned int nbPosition = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 1; T < n; T++)
  {
    ///
    /// FILL u1 and u2 vectors in position T-1
    ///
    theStart = 0;
    while(theStart < (p - 1) && Q[theStart][T - 1] > Q[theStart + 1][T - 1])
      {theStart = theStart + 1;}
    u1[T-1] = theStart;

    theEnd = p - 1;
    while(theEnd > 0 && Q[theEnd][T - 1] > Q[theEnd - 1][T - 1])
      {theEnd = theEnd - 1;}
    u2[T-1] = theEnd;


    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*T -> restricted on each column
      /////
      temp_indState = 0;
      temp_chpt = 0;

      temp_Q = Q[0][0] + cost.slopeCost(states[0], states[v], zero, T, S12P[0][0], S12P[0][T], S12P[1][0], S12P[1][T], S12P[2][0], S12P[2][T]);

      for(unsigned int t = 0; t < T; t++)
      {
        /////
        ///// FIND the minimum of the cost in start state
        /////
        if(t < (T-1))
        {
          theV = cost.vhat(states[v], t, T, S12P[0][t], S12P[0][T], S12P[2][t], S12P[2][T]);
          indexTheV = cost.closestStateIndex(theV, states, p);
        }
        else
        {
          indexTheV = u1[T-1]; // if t = T-1 cost.slopeCost does not depend on u
        }

        ///
        /// explore values between min(u1[t],indexTheV) and max(u2[t],indexTheV)
        ///
        for(unsigned int u = std::min(u1[t],indexTheV); u < std::max(u2[t],indexTheV) + 1; u++) /////explore column of states
        {
          nbPosition = nbPosition + 1; //we explore +1 position
          temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
          if(temp_Q > temp_cost)
          {
            temp_Q = temp_cost;
            temp_indState = u;
            temp_chpt = t;
          }
        }
      }

      /////
      ///// Write response
      /////
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
      Q[v][T] = temp_Q;
    }
  }

  pruning = 2.0*nbPosition/(1.0*p*p*n*(n-1)); //nbPosition seen / nbPosition total

  delete [] u1;
  u1 = NULL;
  delete [] u2;
  u2 = NULL;
}



//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//

void OmegaSN::algoPruning(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;

  unsigned int Tp1;
  unsigned int nm1 = n-1;

  /// PREPROCESSING
  S12P = preprocessing(data);

  double* S1decay = new double[n];
  S1decay[0] = 0;
  for(unsigned int i = 1; i < n; i++){S1decay[i] = S1decay[i-1] + data[i-1];} //cumsum from 0, y_1,...

  double* MAX_Y = new double[n]; //new type of max
  double* MIN_Y = new double[n]; //new type of min
  for(unsigned int i = 0; i < n; i++)
  {
    MIN_Y[i] = 2.0*data[i];
    MAX_Y[i] = 2.0*data[i];
  }

  for(unsigned int i = 0; i < n-1; i++)
  {
    for(unsigned int j = i + 1; j < n-1; j++)
    {
      MIN_Y[i] = std::min(2.0*(S1decay[j+1] - S1decay[i]), MIN_Y[i]);
      MAX_Y[i] = std::max(2.0*(S1decay[j+1] - S1decay[i]), MAX_Y[i]);
    }
  }

  std::list< unsigned int>* t_pos = new std::list< unsigned int>[p];
  std::list< unsigned int>* u_pos = new std::list< unsigned int>[p];
  std::list<unsigned int>::iterator t_it;
  std::list<unsigned int>::iterator u_it;

  Costs cost;
  ///
  /// FILL FIRST COLUMN in Q
  ///
  for(unsigned int i = 0; i < p; i++)
  {
    Q[i][0] = (data[0] - states[i])*(data[0] - states[i]);
    lastChpt[i][0] = 0;
    lastIndState[i][0] = 0;
  }

  ///
  /// ALGO
  ///
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;
  ///variables for pruning
  double delta;
  double DELTA;
  double K;
  int nbPosition = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 1; T < n; T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// Add last column to explore
      /////
      for(unsigned int w = 0; w < p; w++)
      {
        u_pos[v].push_back(w);
        t_pos[v].push_back(T-1);
      }

      /// FIRST ELEMENT
      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();

      temp_Q = Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S12P[0][*t_it], S12P[0][T], S12P[1][*t_it], S12P[1][T], S12P[2][*t_it], S12P[2][T]);
      temp_indState = *u_it;
      temp_chpt = *t_it;

      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();
      while(t_it != t_pos[v].end())
      {
        nbPosition = nbPosition + 1;
        temp_cost = Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S12P[0][*t_it], S12P[0][T], S12P[1][*t_it], S12P[1][T], S12P[2][*t_it], S12P[2][T]);
        if(temp_Q > temp_cost)
        {
          temp_Q = temp_cost;
          temp_indState = *u_it;
          temp_chpt = *t_it;
        }
        ++u_it;
        ++t_it;
      }

      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;

      ////////////////////////
      ///// PRUNING STEP /////
      ////////////////////////
      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();
      while(t_it != t_pos[v].end())
      {
        if(2 < 1)
        {
          u_it = u_pos[v].erase(u_it);
          t_it = t_pos[v].erase(t_it);
        }
        else{++u_it; ++t_it;}
      }

    }
  }

  pruning = 2.0*nbPosition/(1.0*p*p*n*(n-1)); //nbPosition seen / nbPosition total

  delete [] S1decay;
  S1decay = NULL;
  delete [] MAX_Y;
  MAX_Y = NULL;
  delete [] MIN_Y;
  MIN_Y = NULL;
  delete [] t_pos;
  t_pos = NULL;
  delete [] u_pos;
  u_pos = NULL;
}


//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//

void OmegaSN::backtracking(unsigned int n)
{
  unsigned int p = nbStates;
  double temp_Q = Q[0][n-1];
  unsigned int temp_indState = 0;

  for(unsigned int v = 1; v < p; v++)
  {
    if(Q[v][n-1] < temp_Q)
    {
      temp_Q = Q[v][n-1];
      temp_indState = v;
    }
  }

  globalCost = Q[temp_indState][n-1]; ///fill the globalCost OmegaSN variable
  unsigned int temp_chpt = n-1;

  while(temp_chpt > 0)
  {
    changepoints.push_back(temp_chpt);
    parameters.push_back(states[temp_indState]);

    temp_chpt = lastChpt[temp_indState][temp_chpt];
    temp_indState = lastIndState[temp_indState][changepoints[changepoints.size()-1]];
  }

  changepoints.push_back(0);
  parameters.push_back(states[temp_indState]);

  ///reverse the vector
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(parameters.begin(), parameters.end());
}

//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//channel pruning

void OmegaSN::algoISOTONIC(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;
  unsigned int zero = 0;

  /// PREPROCESSING
  S12P = preprocessing(data);

  Costs cost;
  ///
  /// FILL FIRST COLUMN in Q
  ///
  for(unsigned int i = 0; i < p; i++)
  {
    Q[i][0] = (data[0] - states[i])*(data[0] - states[i]);
    lastIndState[i][0] = 0;
    lastChpt[i][0] = 0;
  }

  ///
  /// ALGO
  ///
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///
  /// CHANNEL INFORMATION
  /// u1 / u2 = "min / max" in each column of Q
  ///
  unsigned int* u1 = new unsigned int[n];
  unsigned int* u2 = new unsigned int[n];
  unsigned int theStart;
  unsigned int theEnd;
  double theV = 0;
  unsigned int indexTheV = 0;
  unsigned int nbPosition = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 1; T < n; T++)
  {
    ///
    /// FILL u1 and u2 vectors
    ///
    theStart = 0;
    while(theStart < (p - 1) && Q[theStart][T - 1] > Q[theStart + 1][T - 1])
      {theStart = theStart + 1;}
    u1[T-1] = theStart;

    theEnd = p - 1;
    while(theEnd > 0 && Q[theEnd][T - 1] > Q[theEnd - 1][T - 1])
      {theEnd = theEnd - 1;}
    u2[T-1] = theEnd;


    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*T -> restricted on each column
      /////
      temp_Q = Q[0][0] + cost.slopeCost(states[0], states[v], zero, T, S12P[0][0], S12P[0][T], S12P[1][0], S12P[1][T], S12P[2][0], S12P[2][T]);
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 0; t < T; t++)
      {
        /////
        ///// FIND the minimum of the cost in start state
        /////
        if(t < T-1)
        {
          theV = cost.vhat(states[v], t, T, S12P[0][t], S12P[0][T], S12P[2][t], S12P[2][T]);
          indexTheV = cost.closestStateIndex(theV, states, p);
        }
        else
        {
          indexTheV = u1[T-1];
        }

        ///
        /// explore values between std::min(std::min(u1[t],indexTheV), v); u < std::min(std::max(u2[t],indexTheV), v) + 1
        ///
        for(unsigned int u = std::min(std::min(u1[t],indexTheV), v); u < std::min(std::max(u2[t],indexTheV), v) + 1; u++) /////explore column of states
        {
          nbPosition = nbPosition + 1;
          temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
          if(temp_Q > temp_cost)
          {
            temp_Q = temp_cost;
            temp_indState = u;
            temp_chpt = t;
          }
        }
      }

      /////
      ///// Write response
      /////
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
      Q[v][T] = temp_Q;
    }
  }

  pruning = 2.0*nbPosition/(1.0*p*p*n*(n-1)); //nbPosition seen / nbPosition total

  delete [] u1;
  u1 = NULL;
  delete [] u2;
  u2 = NULL;
}




//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
// NO PRUNING

void OmegaSN::algoUNIMODAL(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;

  /// PREPROCESSING
  S12P = preprocessing(data);

  int** SLOPE = new int*[p];
  for(unsigned int i = 0; i < p; i++){SLOPE[i] = new int[n];}

  Costs cost;
  ///
  /// FILL FIRST COLUMN in Q
  ///
  for(unsigned int i = 0; i < p; i++)
  {
    Q[i][0] = (data[0] - states[i])*(data[0] - states[i]);
    lastChpt[i][0] = 0;
    lastIndState[i][0] = 0;
    SLOPE[i][0] = 0;
  }

  ///
  /// ALGO
  ///
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///


  for(unsigned int T = 1; T < n; T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*T
      /////
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 0; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore column of states
        {
          if(!(u < v && SLOPE[u][t] == -1))
          {
            temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
            if(temp_Q > temp_cost)
            {
              temp_Q = temp_cost;
              temp_indState = u;
              temp_chpt = t;
            }
          }
        }

      }
      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;

      if(temp_indState > v){SLOPE[v][T] = -1;}
      if(temp_indState < v){SLOPE[v][T] = 1;}
      if(temp_indState == v){if(SLOPE[temp_indState][temp_chpt] == -1){SLOPE[v][T] = -1;}}
    }
  }

  pruning = 1;

  for(unsigned int i = 0; i < p; i++){delete(SLOPE[i]);}
  delete [] SLOPE;
  SLOPE = NULL;

}



//####### algoOUTLIER #######////####### algoOUTLIER #######////####### algoOUTLIER #######//
//####### algoOUTLIER #######////####### algoOUTLIER #######////####### algoOUTLIER #######//
//####### algoOUTLIER #######////####### algoOUTLIER #######////####### algoOUTLIER #######//
//####### algoOUTLIER #######////####### algoOUTLIER #######////####### algoOUTLIER #######//
//####### algoOUTLIER #######////####### algoOUTLIER #######////####### algoOUTLIER #######//
// NO PRUNING

void OmegaSN::algoSMOOTHING(std::vector< double >& data, double minAngle)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;

  /// PREPROCESSING
  S12P = preprocessing(data);

  Costs cost;
  ///
  /// FILL FIRST COLUMN in Q
  ///
  for(unsigned int i = 0; i < p; i++)
  {
    Q[i][0] = (data[0] - states[i])*(data[0] - states[i]);
    lastChpt[i][0] = 0;
    lastIndState[i][0] = i;
  }

  ///
  /// ALGO
  ///
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  //std::cout << "minAngle " << minAngle << std::endl;

  for(unsigned int T = 1; T < n; T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*T
      /////
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 0; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore column of states
        {
          //std::cout << "-- " << lastChpt[u][t] << " -- ";
          if(cost.angleTest(lastChpt[u][t], t, T, states[lastIndState[u][t]], states[u], states[v], minAngle))
          {
            temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
            if(temp_Q > temp_cost)
            {
              temp_Q = temp_cost;
              temp_indState = u;
              temp_chpt = t;
            }
          }
        }

      }
      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;
    }
  }
  pruning = 1;

}
