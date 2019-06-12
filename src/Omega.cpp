#include "Omega.h"
#include "Costs.h"

#include <algorithm> // std::reverse // std::min // std::max
#include<iostream>
#include <stdlib.h>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Omega::Omega(std::vector< double >& values, double beta, unsigned int n)
{
  nbStates = values.size();
  states = new double[nbStates];
  for(unsigned int i = 0; i < nbStates; i++){states[i] = values[i];}

  ///
  /// MATRIX INITIALIZATION
  ///
  Q = new double*[nbStates]; ///matrix of costs
  lastChpt = new unsigned int*[nbStates]; ///matrix of best last changepoints
  lastIndState = new unsigned int*[nbStates]; ///matrix of starting states for the best last segment

  for(unsigned int i = 0; i < nbStates; i++){Q[i] = new double[n];}
  for(unsigned int i = 0; i < nbStates; i++){lastChpt[i] = new unsigned int[n];}
  for(unsigned int i = 0; i < nbStates; i++){lastIndState[i] = new unsigned int[n];}

  penalty = beta;
}



Omega::~Omega()
{
  delete(states);
  states = NULL;
  for(unsigned int i = 0; i < nbStates; i++){delete(Q[i]);}
  for(unsigned int i = 0; i < nbStates; i++){delete(lastChpt[i]);}
  for(unsigned int i = 0; i < nbStates; i++){delete(lastIndState[i]);}
  delete [] Q;
  Q = NULL;
  delete [] lastChpt;
  lastChpt = NULL;
  delete [] lastIndState;
  lastIndState = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > Omega::GetChangepoints() const {return(changepoints);}
std::vector< double > Omega::GetParameters() const {return(parameters);}
double Omega::GetGlobalCost() const {return(globalCost);}


//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//
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
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;
  unsigned int zero = 0;

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
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;
    }
  }

  delete(S1);
  S1 = NULL;
  delete(S2);
  S2 = NULL;
  delete(SP);
  SP = NULL;
}

//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
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

  unsigned int zero = 0;
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
    while(theStart < (p - 1) && Q[theStart][T - 1] > Q[theStart + 1][T - 1]){theStart = theStart + 1;}
    u1[T-1] = theStart;
    theEnd = p - 1;
    while(theEnd > 0 && Q[theEnd][T - 1] > Q[theEnd - 1][T - 1]){theEnd = theEnd - 1;}
    u2[T-1] = theEnd;


    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// EXPLORE MATRIX size p*T -> restricted on each column
      /////
      temp_indState = 0;
      temp_chpt = 0;

      temp_Q = Q[0][0] + cost.slopeCost(states[0], states[v], zero, T, S1[0], S1[T], S2[0], S2[T], SP[0], SP[T]) + penalty;


      for(unsigned int t = 0; t < T; t++)
      {
        /////
        ///// FIND the minimum of the cost in start state
        /////
        if(t < T-1){
        theV = cost.vhat(states[v], t, T, S1[t], S1[T], SP[t], SP[T]);
        indexTheV = cost.closestState(theV, states, p);
        }else{indexTheV = u1[T-1];}

        /// std::cout << T << " theV " << theV << "  indexTheV " << indexTheV << "  --u1[t] " << u1[t] << "  u2[t] " << u2[t] << std::endl;
        /// explore values between min(u1[t],indexTheV) and max(u2[t],indexTheV)
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

      /////
      ///// Write response
      /////
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
      Q[v][T] = temp_Q;
    }
  }

  delete(u1);
  u1 = NULL;
  delete(u2);
  u2 = NULL;
  delete(S1);
  S1 = NULL;
  delete(S2);
  S2 = NULL;
  delete(SP);
  SP = NULL;
}



//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//

void Omega::algoPruning(std::vector< double >& data)
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
  SP[0] = data[0];
  for(unsigned int i = 1; i < n; i++){S1[i] = S1[i-1] + data[i];}
  for(unsigned int i = 1; i < n; i++){S2[i] = S2[i-1] + (data[i] * data[i]);}
  for(unsigned int i = 1; i < n; i++){SP[i] = SP[i-1] + (i+1) * data[i];}

  double* MAX_Y = new double[n];
  double* MIN_Y = new double[n];
  MAX_Y[n-1] = 1.0*data[n-2];
  MIN_Y[n-1] = 1.0*data[n-2];
  for(int i = n-2; i > -1; i--){MAX_Y[i] = std::max(1.0*data[i], MAX_Y[i+1]);}
  for(int i = n-2; i > -1; i--){MIN_Y[i] = std::min(1.0*data[i], MIN_Y[i+1]);}


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
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;
  ///variables for pruning
  double delta;
  double DELTA;
  double K;

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

      //std::cout << u_pos[v].size() << " " << t_pos[v].size() << *u_it << *t_it;

      temp_Q = Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S1[*t_it], S1[T], S2[*t_it], S2[T], SP[*t_it], SP[T]) + penalty;
      temp_indState = *u_it;
      temp_chpt = *t_it;

      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();
      while (t_it != t_pos[v].end())
      {
        if(temp_Q > Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S1[*t_it], S1[T], S2[*t_it], S2[T], SP[*t_it], SP[T]) + penalty)
        {
          temp_Q = Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S1[*t_it], S1[T], S2[*t_it], S2[T], SP[*t_it], SP[T]) + penalty;
          temp_indState = *u_it;
          temp_chpt = *t_it;
        }

        ++u_it;
        ++t_it;
      }

      //std::cout << " end "<< *u_it << " "<< *t_it << std::endl;

      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;


      /////
      ///// PRUNING STEP TO BE DONE
      /////
      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();
      while (t_it != t_pos[v].end())
      {
        unsigned int Tp1 = T+1;
        unsigned int nm1 = n-1;
        DELTA = states[*u_it] - states[v];
        if(DELTA >= 0){delta = MAX_Y[T] - states[v];}else{delta = MIN_Y[T] - states[v];}

        K = SP[T] - SP[*t_it] -  (*t_it + 1) * (S1[T] - S1[*t_it]);
        //std::cout << "delta " << delta << " DELTA " << DELTA << " K " << K << std::endl;

        if((Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S1[*t_it], S1[T], S2[*t_it], S2[T], SP[*t_it], SP[T]) > temp_Q) && cost.pruningTest(*t_it, T, Tp1, delta, DELTA, K, states[v]) && cost.pruningTest(*t_it, T, nm1, delta, DELTA, K, states[v]))
          {u_it = u_pos[v].erase(u_it); t_it = t_pos[v].erase(t_it);}else{++u_it; ++t_it;}

      }

    }
  }


  delete(S1);
  S1 = NULL;
  delete(S2);
  S2 = NULL;
  delete(SP);
  SP = NULL;
  delete(MAX_Y);
  MAX_Y = NULL;
  delete(MIN_Y);
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
