// MIT License
// Copyright (c) 2019 Marco Pascucci

#include <assert.h>
#include <iostream>
#include <vector>

using namespace std;

template <typename T> void printv(vector<T> vector) {
  cout << "[ ";
  for (size_t i = 0; i < vector.size(); i++) {
    cout << vector[i];
    if (i != vector.size() - 1)
      cout << ", ";
  }
  cout << " ]" << endl;
}

template <class T> vector<T> range(int start, int end) {
  if (start > end) {
    throw std::invalid_argument("start > end");
  }
  vector<T> out(end - start);
  for (size_t i = 0; i < end - start; i++) {
    out[i] = i + start;
  }
  return out;
}

template <class T> vector<T> range(int end) { return range<T>(0, end); }

template <class T> void range(vector<T> &v) {
  for (typename vector<T>::size_type i = 0; i < v.size(); i++) {
    v[i] = i;
  }
}

// template<typename T1, typename T2>
// double dot(vector<T1> &a, vector<T2> &b, size_t start, size_t end) {
//     if(a.size() != b.size()){
//         throw std::invalid_argument( "a and b must have the same dimension"
//         );
//     }
//     if(start < 0 || end > a.size() ){
//         throw std::invalid_argument( "invalid start/end range" );
//     }
//     double sum = 0;
//     for (size_t i=start; i< end; i++) {
//         sum = sum + a[i]*b[i];
//     }
//     return sum;
// }
// template<typename T1, typename T2>
// double dot(vector<T1> &a, vector<T2> &b) {
//     return dot(a, b, 0, a.size());
// }
//  //specialization to return int
// int dot(vector<int> &a, vector<int> &b, size_t start, size_t end);
// int dot(vector<int> &a, vector<int> &b);

// template<class T>
// T accumulate(vector<T> &v, size_t start, size_t end) {
//     if(start < 0 || end > v.size() ){
//         throw std::invalid_argument( "invalid start/end range" );
//     }
//     T sum = 0;
//     for (typename vector<T>::size_type i=start; i<end; i++) {
//         sum = sum + v[i];
//     }
//     return sum;
// }
// template<class T>
// T accumulate(vector<T> &v) {
//     return accumulate(v, 0, v.size());
// }

template <typename T1, typename T2>
void lin_reg_no_const(vector<T1> &x, vector<T2> &y, double *coeff, size_t start,
                      size_t end) {
  *coeff =
      1.0 /
      (inner_product(x.begin() + start, x.begin() + end, x.begin() + start,
                     0.0)) *
      inner_product(x.begin() + start, x.begin() + end, y.begin() + start, 0.0);
}

template <typename T1, typename T2>
void lin_reg(vector<T1> &x, vector<T2> &y, double *coeff, double *intercept,
             size_t start, size_t end) {
  /* calculate LSLR on x,y from position start to end (excluded). */

  assert(y.size() == x.size());
  assert(start >= 0 && end <= x.size());

  double m00, m11, m01, det;
  double i00, i11, i01;

  m00 = end - start;
  m11 =
      inner_product(x.begin() + start, x.begin() + end, x.begin() + start, 0.0);
  m01 = accumulate(x.begin() + start, x.begin() + end, 0.0);

  det = (m00 * m11) - m01 * m01;

  i00 = m11 / det;
  i11 = m00 / det;
  i01 = -m01 / det;

  double v1 = accumulate(y.begin() + start, y.begin() + end, 0.0);
  double v2 =
      inner_product(x.begin() + start, x.begin() + end, y.begin() + start, 0.0);

  *intercept = i00 * v1 + i01 * v2;
  *coeff = i01 * v1 + i11 * v2;
}

template <typename T1, typename T2>
void lin_reg(vector<T1> &x, vector<T2> &y, double *coeff, double *intercept) {
  return lin_reg(x, y, coeff, intercept, 0, x.size());
}
template <typename T1, typename T2>
void lin_reg_no_const(vector<T1> &x, vector<T2> &y, double *coeff) {
  return lin_reg_no_const(x, y, coeff, 0, x.size());
}

template <typename T>
vector<double> line(const vector<T> &x, const double coeff, size_t start,
                    size_t end, const double inter = 0) {
  /* calculate a linear segment, given a vector x, and the line coefficients
   * from start to end (excluded)
   */

  if (start < 0 || end > x.size()) {
    throw std::invalid_argument("invalid start/end range");
  }
  vector<double> y(end - start);
  for (size_t i = start; i < end; i++) {
    y[i] = coeff * x[i] + inter;
  }
  return y;
}
template <typename T>
vector<double> line(const vector<T> &x, const double coeff,
                    const double inter = 0) {
  return line(x, coeff, 0, x.size(), inter);
}
