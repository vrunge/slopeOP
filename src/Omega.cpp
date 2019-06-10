#include "Omega.h"
#include "Costs.h"

#include <algorithm>    // std::reverse // std::min // std::max
#include<iostream>
#include <stdlib.h>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Omega::Omega(std::vector< double >& values, double beta)
{
  nbStates = values.size();
  states = new double[nbStates];
  for(unsigned int i = 0; i < nbStates; i++){states[i] = values[i];}

  Q = new double*[nbStates]; ///matrix of cost
  lastChpt = new unsigned int*[nbStates]; ///matrix of cost
  lastIndState = new unsigned int*[nbStates]; ///matrix of cost

  penalty = beta;
}

Omega::~Omega()
{
  delete [] states;
  states = NULL;

  for(unsigned int i = 0; i < nbStates; i++){delete [] Q[i]; Q[i] = NULL;}
  for(unsigned int i = 0; i < nbStates; i++){delete [] lastChpt[i]; lastChpt[i] = NULL;}
  for(unsigned int i = 0; i < nbStates; i++){delete [] lastIndState[i]; lastIndState[i] = NULL;}
}


//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > Omega::GetChangepoints() const {return(changepoints);}
std::vector< double > Omega::GetParameters() const {return(parameters);}
double Omega::GetGlobalCost() const {return(globalCost);}


//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//

void Omega::algo(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;

  ///
  /// PREPROCESSING
  ///
  double* S1 = new double[n];
  double* S2 = new double[n];
  double* SP = new double[n];
  S1[0] = data[0];
  S2[0] = data[0] * data[0];
  SP[0] = data[0] ;
  for(unsigned int i = 1; i < n; i++){S1[i] = S1[i-1] + data[i];}
  for(unsigned int i = 1; i < n; i++){S2[i] = S2[i-1] + (data[i] * data[i]);}
  for(unsigned int i = 1; i < n; i++){SP[i] = SP[i-1] + (i+1) * data[i];}

  ///
  /// MATRIX INITIALIZATION
  ///
  for(unsigned int i = 0; i < p; i++){Q[i] = new double[n];}
  for(unsigned int i = 0; i < p; i++){lastChpt[i] = new unsigned int[n];}
  for(unsigned int i = 0; i < p; i++){lastIndState[i] = new unsigned int[n];}

  Costs cost;
  ///
  /// FIRST STEP
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
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;
  unsigned int zero = 0;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  /// states u to v /// time position t to T

  for(unsigned int T = 1; T < n; T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*j
      /////
      temp_Q = Q[0][0] + cost.slopeCost(states[0], states[v], zero, T, S1[0], S1[T], S2[0], S2[T], SP[0], SP[T]) + penalty;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 0; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore colum of states
        {
          if(temp_Q > Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S1[t], S1[T], S2[t], S2[T], SP[t], SP[T]) + penalty)
          {
            temp_Q = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S1[t], S1[T], S2[t], S2[T], SP[t], SP[T]) + penalty;
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

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  delete [] S1;
  delete [] S2;
  delete [] SP;
}


//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//

void Omega::algoChannel(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;

  ///
  /// PREPROCESSING
  ///
  double* S1 = new double[n];
  double* S2 = new double[n];
  double* SP = new double[n];
  S1[0] = data[0];
  S2[0] = data[0] * data[0];
  SP[0] = data[0] ;
  for(unsigned int i = 1; i < n; i++){S1[i] = S1[i-1] + data[i];}
  for(unsigned int i = 1; i < n; i++){S2[i] = S2[i-1] + (data[i] * data[i]);}
  for(unsigned int i = 1; i < n; i++){SP[i] = SP[i-1] + (i+1) * data[i];}

  ///
  /// MATRIX INITIALIZATION
  ///
  for(unsigned int i = 0; i < p; i++){Q[i] = new double[n];}
  for(unsigned int i = 0; i < p; i++){lastChpt[i] = new unsigned int[n];}
  for(unsigned int i = 0; i < p; i++){lastIndState[i] = new unsigned int[n];}

  Costs cost;
  ///
  /// FIRST STEP
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
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  unsigned int* u1 = new unsigned int[n];
  unsigned int* u2 = new unsigned int[n];
  unsigned int theStart;
  unsigned int theEnd;
  double theV = 0;
  unsigned int indexTheV = 0;
  unsigned int zero = 0;

  /// states u to v /// time position t to T

  for(unsigned int T = 1; T < n; T++)
  {
    /// FILL u1 and u2 vectors
    theStart = 0;
    while(theStart < (p - 1) && Q[theStart][T - 1] > Q[theStart + 1][T - 1]){theStart = theStart + 1;}
    u1[T-1] = theStart;
    theEnd = p - 1;
    while(theEnd > 0 && Q[theEnd][T - 1] > Q[theEnd - 1][T - 1]){theEnd = theEnd - 1;}
    u2[T-1] = theEnd;


    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*j
      /////
      temp_indState = 0;
      temp_chpt = 0;

      temp_Q = Q[0][0] + cost.slopeCost(states[0], states[0], zero, T, S1[0], S1[T], S2[0], S2[T], SP[0], SP[T]) + penalty;


      for(unsigned int t = 0; t < T-1; t++)
      {
        theV = cost.vhat(states[v], t, T, S1[t], S1[T], SP[t], SP[T]);
        indexTheV = cost.closestState(theV, states, p);

        //std::cout << T << " theV " << theV << "  indexTheV " << indexTheV << "  --u1[t] " << u1[t] << "  u2[t] " << u2[t] << std::endl;

        for(unsigned int u = std::min(u1[t],indexTheV); u < std::max(u2[t],indexTheV) + 1; u++) /////explore colum of states
        {
          if(temp_Q > Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S1[t], S1[T], S2[t], S2[T], SP[t], SP[T]) + penalty)
          {
            temp_Q = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S1[t], S1[T], S2[t], S2[T], SP[t], SP[T]) + penalty;
            temp_indState = u;
            temp_chpt = t;
          }
        }
      }
      //std::cout << T << "  u1[t] " << u1[T-1] << "  u2[t] " << u2[T-1] << std::endl;

      unsigned int Tm1 = T-1;
      for(unsigned int u = u1[Tm1]; u < u2[Tm1]+1; u++) /////explore colum of states
      {
        if(temp_Q > Q[u][T-1] + cost.slopeCost(states[u], states[v], Tm1, T, S1[Tm1], S1[T], S2[Tm1], S2[T], SP[Tm1], SP[T]) + penalty)
        {
          temp_Q = Q[u][T-1] + cost.slopeCost(states[u], states[v], Tm1, T, S1[Tm1], S1[T], S2[Tm1], S2[T], SP[Tm1], SP[T]) + penalty;
          temp_indState = u;
          temp_chpt = Tm1;
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

  delete [] u1;
  delete [] u2;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  delete [] S1;
  delete [] S2;
  delete [] SP;
}




//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//

void Omega::backtracking(unsigned int n)
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

  globalCost = Q[temp_indState][n-1];
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
