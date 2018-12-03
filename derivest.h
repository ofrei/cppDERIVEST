#ifndef DERIVEST_H__
#define DERIVEST_H__

#include <functional>

enum DerivestStyle {
  DerivestStyle_Central = 0,
  DerivestStyle_Forward = 1,
  DerivestStyle_Backward = 2,
};

// DERIVEST calculates numeric derivative of a function.
// Arguments (input)
//   fun              - function to differentiate
//   x0               - point at which to differentiate fun
//   derivative_order - specifies the derivative order estimated. Must be a positive integer from the set[1, 2, 3, 4].
//   method_order     - specifies the order of the basic method used for the estimation.
//                      for central methods, methods_order must be 2 or 4; otherwise can be 1, 2, 3 or 4
//   style            - specifies the style of the basic method used for the estimation. 'central', 'forward', or 'backwards' difference methods are used.
//   romberg_terms    - Allows the user to specify the generalized Romberg extrapolation method used, or turn it off completely.
//                      Must be a positive integer from the set[0, 1, 2, 3].
//
// Arguments (output)
//   *der             - pointer to where store the result (derivative)
//   *err             - pointer to where store error estimate
//   *finaldelta      - pointer to where store the final overall stepsize chosen by cppDERIVEST
//
// Return value: true if succeeded, false otherwise.
//
// Remark. Recommended params are as follows:
//   style = DerivestStyle_Central
//   method_order = 4
//   romberg_terms = 2
bool derivest(std::function<double(double)> fun, double x0, int derivative_order, int method_order, DerivestStyle style, int romberg_terms, double *der, double* err, double* finaldelta);

// estimate elements of the Hessian matrix (matrix of 2nd partials)
// n gives the dimension
// hess and err are output arrays, must be at least of size n*n
bool hessian(std::function<double(double*)> fun, double* x0, int n, double* hess, double* err);

#endif