# cppDERIVEST
C++ port of [Adaptive Robust Numerical Differentiation](https://se.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation). The algorithm is described [here](http://convexoptimization.com/TOOLS/DERIVEST.pdf). 

# Getting started
Include ``derivest.h`` and compile ``derivest.cc`` with your application.

# Notes

* I've only tested this as ``C++`` code, but with minor changes it should work as plain ``C`` code. If you fix this, please submit pull request.
* Functionality is slightly limited compared to the original (MATLAB) version. Most importantly, ``StepRatio`` parameter is hardcoded to ``2.0000001``

# Demo
```
#include "derivest.h"
#include <math.h>
#include <stdlib.h>
int main() {
  double der=0.0, err = 0.0, finaldelta = 1.0;
  printf("f(x) = exp(x); x0 = 1;\n");
  derivest(exp, 1, 1, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f'(x)    = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - exp(1), err);
  derivest(exp, 1, 2, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f''(x)   = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - exp(1), err);
  derivest(exp, 1, 3, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f'''(x)  = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - exp(1), err);
  derivest(exp, 1, 4, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f''''(x) = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - exp(1), err);
  printf("\nf(x) = x + x*x + x*x*x + x*x*x*x + x*x*x*x*x; x0 = 0.0;\n");
  derivest(poly, 0, 1, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f'(x)    = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - 1, err);
  derivest(poly, 0, 2, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f''(x)   = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - 2, err);
  derivest(poly, 0, 3, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f'''(x)  = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - 6, err);
  derivest(poly, 0, 4, 4, DerivestStyle_Central, 2, &der, &err, NULL); printf("f''''(x) = %.17f, err=%.3e, hat(err)=%.3e\n", der, der - 24, err);
}
```
Produce the following result:
```
f(x) = exp(x); x0 = 1;
f'(x)    = 2.71828182845905353, err=8.438e-15, hat(err)=1.078e-13
f''(x)   = 2.71828182846117983, err=2.135e-12, hat(err)=2.294e-12
f'''(x)  = 2.71828182843829991, err=-2.075e-11, hat(err)=4.355e-11
f''''(x) = 2.71828130980468119, err=-5.187e-07, hat(err)=5.830e-06

f(x) = x + x*x + x*x*x + x*x*x*x + x*x*x*x*x; x0 = 0.0;
f'(x)    = 1.00000000000000000, err=0.000e+00, hat(err)=7.998e-15
f''(x)   = 2.00000000000000355, err=3.553e-15, hat(err)=6.035e-14
f'''(x)  = 6.00000000000983036, err=9.830e-12, hat(err)=6.538e-11
f''''(x) = 23.99999999999797495, err=-2.025e-12, hat(err)=1.920e-11
```
