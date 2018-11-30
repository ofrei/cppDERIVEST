#include <assert.h>

#include <iostream>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_sf_bessel.h>

// 'DerivativeOrder', 'MethodOrder', 'Style', 'RombergTerms', 'MaxStep'

enum DerivestStyle {
  DerivestStyle_Central = 0,
  DerivestStyle_Forward = 1,
  DerivestStyle_Backward = 2,
};

// Matrix with both dimentions no larger than 10. User must know what the actual shape is.
struct SmallMatrix {
public:
  SmallMatrix() { for (int i = 0; i < 100; i++) vals[i] = 0; }
  double vals[100];

  double get(int row, int col) { assert(row < 10 && col < 10); return vals[row + 10 * col]; }
  void set(int row, int col, double val) { assert(row < 10 && col < 10); vals[row + 10 * col] = val; }
  void show() {
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 10; j++) printf("%.4f\t", get(i, j));
      printf("\n");
    }
  }
};

void fdamat(double sr, int parity, int nterms, SmallMatrix* mat) {
  // Compute matrix for fda derivation.
  // parity can be
  //   0 (one sided, all terms included but zeroth order)
  //   1 (only odd terms included)
  //   2 (only even terms included)
  // nterms - number of terms

  // sr is the ratio between successive steps
  
  const static double factorial_case0[] = { 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };  // factorial(1:nterms)
  const static double factorial_case1[] = { 1,    6,     120,      5040,        362880,           39916800, 6227020800, 1307674368000, 355687428096000, 121645100408832000 };  // factorial(1:2:(2*nterms))
  const static double factorial_case2[] = {    2,    24,      720,       40320,         3628800,  479001600, 87178291200, 20922789888000, 6402373705728000, 2432902008176640000 };  // factorial(2:2:(2*nterms))

  const double srinv = 1.0 / sr;

  assert(parity >= 0 & parity <= 2);
  if (parity==0) // single sided rule
    for (int i = 1; i <= nterms; i++)
      for (int j = 1; j <= nterms; j++)
        mat->set(i-1, j-1, 1.0 / factorial_case0[j-1] * pow(srinv, (i - 1)*j));  // c(j).*srinv.^((i-1).*j);
  if (parity==1) // odd order derivative
    for (int i = 1; i <= nterms; i++)
      for (int j = 1; j <= nterms; j++)
        mat->set(i-1, j-1, 1.0 / factorial_case1[j-1] * pow(srinv, (i - 1)*(2*j - 1)));  // c(j).*srinv.^((i-1).*(2*j-1));
  if (parity == 2)  // even order derivative
    for (int i = 1; i <= nterms; i++)
      for (int j = 1; j <= nterms; j++)
        mat->set(i-1, j-1, 1.0 / factorial_case2[j-1] * pow(srinv, (i - 1)*(2*j)));  // c(j).*srinv.^((i-1).*(2*j));
} // fdamat

struct DerivestParams {
public:
  int derivative_order;    // can be 1,2,3,4
  int method_order;        // 2 or 4 for central; 1,2,3 or 4 for forward or backward
  DerivestStyle style;
  int romberg_terms;       // can be 0,1,2,3
  double max_step;
  double step_ratio;

  DerivestParams() : 
      derivative_order(1), 
      method_order(4), 
      style(DerivestStyle_Central),
      romberg_terms(2), 
      max_step(100),
      step_ratio(2.0000001) {
    if ((derivative_order < 1) || (derivative_order > 4)) throw std::invalid_argument("derivative_order must be in 1,2,3,4");

    switch (style) {
    case DerivestStyle_Central: 
      if (method_order != 2 && method_order != 4) throw std::invalid_argument("For style==central method_order must be 2 or 4");
      break;       
    case DerivestStyle_Forward:
    case DerivestStyle_Backward:
      if (method_order < 1 && method_order > 4) throw std::invalid_argument("For style==forward and style==backward method_order must be in 1,2,3,4");
      break;
    default: 
      throw std::invalid_argument("style is misspecified");
    }

    if ((romberg_terms < 0) || (romberg_terms > 3)) throw std::invalid_argument("romberg_terms must be in 0,1,2,3");
    if (max_step < 0) throw std::invalid_argument("max_step must be >= 0");
    if (step_ratio < 0) throw std::invalid_argument("step_ratio must be >= 0");
  }

  void Get() {

  }

};

double derivest(double(*func)(double), double x0, DerivestParams par, double* err, double* finaldelta) {
  const double nominal_step = std::max(x0, 0.02);
  static const int ndel = 26;
  double delta[ndel];
  for (int i = 0; i < ndel; i++) delta[i] = par.max_step*pow(par.step_ratio, -i);

  return 0.0;
}

int main() {
  derivest(exp, 0, DerivestParams(), nullptr, nullptr);
  SmallMatrix mat;
  fdamat(2.0000001, 0, 10, &mat); mat.show(); printf("\n");
  fdamat(2.0000001, 1, 10, &mat); mat.show(); printf("\n");
  fdamat(2.0000001, 2, 10, &mat); mat.show(); printf("\n");

  return 0;
}