// MIT License
// Copyright (c) 2019 Marco Pascucci

#include <vector>
#include <numeric>
#include <iostream>
#include <assert.h>
#include <cmath>
#include "linalg.hpp"
#include <algorithm>


template<typename Tx, typename Ty>
struct Point {
    Tx x;
    Ty y;
    Point();
    Point(Tx x, Ty y) : x(x), y(y) {};
    Point(Point<Tx,Ty> &p) : x(p.x), y(p.y) {};
    void operator=(Point<Tx,Ty> &p) {
        x = p.x;
        y = p.y;
    }
};


template<typename T1, typename T2>
double sum_of_squares(vector<T1> &y, vector<T2> &y0, size_t start, size_t end) {
    /* @brief Calculate sum of square errors between y and y0 */
    if(y.size() != y0.size()){
        throw std::invalid_argument( "a and b must have the same dimension" );
    }
    if(start < 0 || end > y.size() ){
        throw std::invalid_argument( "invalid start/end range" );
    }
    double sum_of_squares = 0;
    for (size_t i=start; i<end; i++) {
        sum_of_squares = sum_of_squares + pow((double)(y[i]-y0[i]),2);
    }
    return sum_of_squares;
}
template<typename T1, typename T2>
double sum_of_squares(vector<T1> &y, vector<T2> &y0) {
    return sum_of_squares(y,y0,0,y.size());
}


template<typename T1, typename T2>
double cost_linear(vector<T1> &x, vector<T2> &y,
                    double coeff, double inter,
                    size_t start, size_t end,
                    bool average=false) {
    /* calculate the sum of square errors between y and x*coeff+inter */
    assert(y.size() == x.size());
    assert(start >= 0 && end <= x.size());

    double cost = 0;
    for (uint i=start; i<end; i++) {
        cost = cost + pow((y[i] - (x[i]*coeff+inter)),2);
    }
    if (average) cost = cost/x.size();
    return cost;
}
template<typename T1, typename T2>
double cost_linear(vector<T1> &x, vector<T2> &y,
                    double coeff, double inter,
                    bool average=false) {
    return cost_linear(x, y, coeff, inter, 0, x.size(), average);
}


template<typename T1, typename T2>
double cost_linear_point(vector<T1> &x, vector<T2> &y,
                    Point<T1,T2> p,
                    double coeff,
                    size_t start, size_t end,
                    bool average=false) {
    /* calculate the sum of square errors between y and x*coeff+inter */
    assert(y.size() == x.size());
    assert(start >= 0 && end <= x.size());

    double cost = 0;
    for (uint i=start; i<end; i++) {
        cost = cost + pow((y[i] - ((x[i]-p.x)*coeff + p.y)),2);
    }
    if (average) cost = cost/x.size();
    return cost;
}
template<typename T1, typename T2>
double cost_linear_point(vector<T1> &x, vector<T2> &y,
                    Point<T1,T2> p,
                    double coeff,
                    bool average=false) {
    return cost_linear_point(x, y, p, coeff, 0, x.size(), average);
}


vector<uint> backtrack(vector<uint> cp_raw){
    /* @brief  backtrack the changepoint vector */
    uint x = cp_raw[cp_raw.size()-1];
    vector<uint> cp;
    while (x > 0) {
        cp.push_back(x);
        x = cp_raw[x];
    }
    reverse(cp.begin(),cp.end());
    return cp;
}


template<typename Tx, typename Ty>
vector<uint> pelt(vector<Tx> &x, vector<Ty> &y, double beta) {
    /* detect changepoints with PELT algorithm */
    assert(x.size() == y.size());
    int n = y.size();
    vector<double> q(n);
    vector<uint> cp(n,-1);
    vector<uint> P = {0};
    vector<uint> temp;

    double coeff,inter,cost,sum_of_squares,this_cost,q_temp,cp_temp;
    vector<double> q_new;

    cp[0] = 0;

    for (size_t t=1; t<n; ++t) {
        // possible changepoint @ t

        lin_reg(x,y,&coeff,&inter,0,t+1);
        q_temp = cost_linear(x,y,coeff,inter,0,t+1);
        cp_temp = 0;
        q_new = vector<double>(t,0);
        q_new[0] = q_temp;

        // return vector<int>(1,1);

        for (auto tau: P) {
            if (t-tau == 1) continue;
            lin_reg(x,y,&coeff,&inter,tau+1,t+1);
            this_cost = cost_linear(x,y,coeff,inter,tau+1,t+1);
            q_new[tau] = q[tau] + this_cost; // q_new
            if (q_new[tau] + beta < q_temp) {
                q_temp = q_new[tau] + beta;
                cp_temp = tau;
            }
        }
        q[t] = q_temp;
        cp[t] = cp_temp;

        temp.clear();

        for (size_t i=0; i<P.size(); i++) {
            if (q_new[P[i]] < q_temp) {
                // printf("erased: %d", i);
                temp.push_back(P[i]);
            }
        }
        temp.push_back(t);
        P = temp;
        // printv(P);

    }
    return cp;
}


template<typename Tx, typename Ty>
vector<uint> peltcc(vector<Tx> &x, vector<Ty> &y, double beta) {
    /* detect changepoints with PELT algorithm */
    assert(x.size() == y.size());
    int n = y.size();
    vector<double> q(n);
    vector<uint> cp(n,-1);
    vector<uint> P = {0};
    vector<uint> temp;
    vector<Tx> x_trans(y.size());
    vector<Ty> y_trans(y.size());
    vector<Ty> y_temp(y.size());

    double coeff,coeff_temp,inter,cost,sum_of_squares,this_cost,q_temp,cp_temp;
    Point<Tx,Ty> pt_temp(0,0), pt(0,0);
    vector<double> q_new;

    cp[0] = 0;

    for (size_t t=1; t<n; ++t) {
        // possible changepoint @ t

        lin_reg(x,y,&coeff,&inter,0,t+1);
        q_temp = cost_linear(x,y,coeff,inter,0,t+1);
        cp_temp = 0;
        q_new = vector<double>(t,0);
        q_new[0] = q_temp;
        coeff_temp = coeff;
        pt.x = x[0];
        pt.y = coeff*x[0] + inter;
        pt_temp = pt;

        // return vector<int>(1,1);

        for (auto tau: P) {
            if (t-tau == 1) continue;
            pt.x = x[tau];
            pt.y =  y_temp[tau];
            transform(x.begin(), x.end(), x_trans.begin(), bind2nd(std::plus<int>(), -pt.x));
            transform(y.begin(), y.end(), y_trans.begin(), bind2nd(std::plus<double>(), -pt.y));
            lin_reg_no_const(x_trans, y_trans, &coeff, tau+1, t+1);
            this_cost = cost_linear_point(x,y,pt,coeff,tau+1,t+1);
            q_new[tau] = q[tau] + this_cost; // q_new
            if (q_new[tau] + beta < q_temp) {
                q_temp = q_new[tau] + beta;
                cp_temp = tau;
                coeff_temp = coeff;
                pt_temp = pt;
            }
        }
        q[t] = q_temp;
        cp[t] = cp_temp;
        y_temp[t] = coeff_temp*(t - pt_temp.x) + pt_temp.y;

        temp.clear();

        for (size_t i=0; i<P.size(); i++) {
            if (q_new[P[i]] < q_temp) {
                // printf("erased: %d", i);
                temp.push_back(P[i]);
            }
        }
        temp.push_back(t);
        P = temp;
        // printv(P);

    }
    return cp;
}
