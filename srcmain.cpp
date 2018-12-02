#include <assert.h>

#include <iostream>
#include <algorithm>
#include <cmath>


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

double* get_fdarule_test(int i) {
  if (i==0) { static double data[] = { 1.0, 1.2, 1.5 }; return data; }
}

void get_fdarule(int derivative_order, int method_order, DerivestStyle style, double* fdarule, int *n) {
  // WARNING (!) all constants in get_fdarule are hard-coded to a specific step_ratio.
  // const double step_ratio = 2.0000001;  

  *n = 0;
  if (style == DerivestStyle_Central) {
    // for central rules, we will reduce the load by an even or odd transformation as appropriate.
    if (method_order == 2 && derivative_order == 1) { fdarule[(*n)++] = 1; } // the odd transformation did all the work
    else if (method_order == 2 && derivative_order == 2) { fdarule[(*n)++] = 2; } // the even transformation did all the work
    else if (method_order == 2 && derivative_order == 3) { fdarule[(*n)++] = 7.9999997333333637; fdarule[(*n)++] = -16.000000266666699; } // the odd transformation did most of the work, but we need to kill off the linear term
    else if (method_order == 2 && derivative_order == 4) { fdarule[(*n)++] = 31.999998933333455; fdarule[(*n)++] = -128.00000853333367; } // the even transformation did most of the work, but we need to kill off the quadratic term
    else if (method_order == 4 && derivative_order == 1) { fdarule[(*n)++] = -0.33333328888889396; fdarule[(*n)++] = 2.6666667111111164; } // the odd transformation did most of the work, but we need to kill off the cubic term
    else if (method_order == 4 && derivative_order == 2) { fdarule[(*n)++] = -0.66666657777778793; fdarule[(*n)++] = 10.666667377777806; } // the even transformation did most of the work, but we need to kill off the quartic term
    else if (method_order == 4 && derivative_order == 3) { fdarule[(*n)++] = -2.6666662222222755; fdarule[(*n)++] = 90.666673155556012; fdarule[(*n)++] = -170.66668942222367; } // the odd transformation did much of the work, but we need to kill off the linear & quintic terms
    else if (method_order == 4 && derivative_order == 4) { fdarule[(*n)++] = -10.66666488888913; fdarule[(*n)++] = 725.33342151111754; fdarule[(*n)++] = -2730.6673038222893; } // the even transformation did much of the work, but we need to kill off the quadratic and 6th order terms
    else std::invalid_argument("Given combination of derivative_order, method_order and style is not supported");
    return;
  }
  if (style == DerivestStyle_Forward || style == DerivestStyle_Backward) { // These two cases are identical, except at the very end, where a sign will be introduced.
    // No odd/even trans, but we already dropped off the constant term
    // method_order methods drop off the lower order terms, plus terms directly above derivative_order
    if ((method_order = 1) && (derivative_order == 1)) { fdarule[(*n)++] = 1.00000000000000000000; }
    else if ((method_order = 1) && (derivative_order == 2)) { fdarule[(*n)++] = 3.99999980000002120000; fdarule[(*n)++] = -8.00000000000002130000; }
    else if ((method_order = 1) && (derivative_order == 3)) { fdarule[(*n)++] = 15.99999866666682300000; fdarule[(*n)++] = -96.00000000000036900000; fdarule[(*n)++] = 128.00000853333381000000; }
    else if ((method_order = 1) && (derivative_order == 4)) { fdarule[(*n)++] = 73.14284948027338400000; fdarule[(*n)++] = -1024.00001706667810000000; fdarule[(*n)++] = 4096.00047786673080000000; fdarule[(*n)++] = -4681.14377108039120000000; }
    else if ((method_order = 2) && (derivative_order == 1)) { fdarule[(*n)++] = -0.99999990000001060000; fdarule[(*n)++] = 4.00000000000001070000; }
    else if ((method_order = 2) && (derivative_order == 2)) { fdarule[(*n)++] = -3.99999940000007470000; fdarule[(*n)++] = 39.99999920000021100000; fdarule[(*n)++] = -64.00000320000023600000; }
    else if ((method_order = 2) && (derivative_order == 3)) { fdarule[(*n)++] = -15.99999706666723800000; fdarule[(*n)++] = 351.99999466667157000000; fdarule[(*n)++] = -1664.00014933335680000000; fdarule[(*n)++] = 2048.00034133337610000000; }
    else if ((method_order = 2) && (derivative_order == 4)) { fdarule[(*n)++] = -73.14284216598025500000; fdarule[(*n)++] = 3364.57144912082190000000; fdarule[(*n)++] = -36864.00552959775100000000; fdarule[(*n)++] = 135753.17825549439000000000; fdarule[(*n)++] = -149796.62314405275000000000; }
    else if ((method_order = 3) && (derivative_order == 1)) { fdarule[(*n)++] = 0.33333325555556748000; fdarule[(*n)++] = -3.99999960000004600000; fdarule[(*n)++] = 10.66666684444448300000; }
    else if ((method_order = 3) && (derivative_order == 2)) { fdarule[(*n)++] = 1.33333295555564010000; fdarule[(*n)++] = -34.66666284444525600000; fdarule[(*n)++] = 234.66667484444733000000; fdarule[(*n)++] = -341.33337315556037000000; }
    else if ((method_order = 3) && (derivative_order == 3)) { fdarule[(*n)++] = 5.33333164444366050000; fdarule[(*n)++] = -287.99997119996817000000; fdarule[(*n)++] = 4309.33369031081380000000; fdarule[(*n)++] = -18432.00368639926800000000; fdarule[(*n)++] = 21845.33952284409700000000; }
    else if ((method_order = 3) && (derivative_order == 4)) { fdarule[(*n)++] = 24.38094413806237200000; fdarule[(*n)++] = -2681.90455978352110000000; fdarule[(*n)++] = 84065.53641580433800000000; fdarule[(*n)++] = -831683.30242319033000000000; fdarule[(*n)++] = 2946000.48652337310000000000; fdarule[(*n)++] = -3195661.82635462100000000000; }
    else if ((method_order = 4) && (derivative_order == 1)) { fdarule[(*n)++] = -0.04761902834467632300; fdarule[(*n)++] = 1.33333302222228410000; fdarule[(*n)++] = -10.66666577777794500000; fdarule[(*n)++] = 24.38095348390041300000; }
    else if ((method_order = 4) && (derivative_order == 2)) { fdarule[(*n)++] = -0.19047610385484859000; fdarule[(*n)++] = 11.04761640453383700000; fdarule[(*n)++] = -191.99999039998590000000; fdarule[(*n)++] = 1121.52392852603750000000; fdarule[(*n)++] = -1560.38125702673600000000; }
    else if ((method_order = 4) && (derivative_order == 3)) { fdarule[(*n)++] = -0.76190439003349786000; fdarule[(*n)++] = 89.90474156592972600000; fdarule[(*n)++] = -3248.76192546417310000000; fdarule[(*n)++] = 42032.77030951646200000000; fdarule[(*n)++] = -171641.96048256109000000000; fdarule[(*n)++] = 199728.84417419793000000000; }
    else if ((method_order = 4) && (derivative_order == 4)) { fdarule[(*n)++] = -3.48299142217729240000; fdarule[(*n)++] = 828.95221520044026000000; fdarule[(*n)++] = -61049.90956226736300000000; fdarule[(*n)++] = 1656010.53460411840000000000; fdarule[(*n)++] = -15628783.09987751400000000000; fdarule[(*n)++] = 54326255.84047879300000000000; fdarule[(*n)++] = -58434969.54424954200000000000; }
    else std::invalid_argument("Given combination of derivative_order, method_order and style is not supported");

    /* The above code is generated in matlab, using subfunctions from derivest.m from DERIVESTsuite.
    for mo = 1:4, for do = 1:4
      v = zeros(1,do + mo - 1); v(do) = 1; fdarule = v/fdamat(2.0000001,0,do+mo-1);
      fprintf('else if ((method_order = %i) && (derivative_order == %i)) { %s }\n', mo,do,sprintf('fdarule[(*n)++] = %.20f; ', fdarule));
    end; end 
    */

    if (style == DerivestStyle_Backward) for (int i = 0; i < (*n); i++) fdarule[i] *= -1; // correct sign for the 'backward' rule
    return;
  }

  throw std::invalid_argument("style is misspecified");
}

struct DerivestParams {
public:
  int derivative_order;    // can be 1,2,3,4
  int method_order;        // 2 or 4 for central; 1,2,3 or 4 for forward or backward
  DerivestStyle style;
  int romberg_terms;       // can be 0,1,2,3
  double max_step;

  DerivestParams() : 
      derivative_order(1), 
      method_order(4), 
      style(DerivestStyle_Central),
      romberg_terms(2), 
      max_step(100) {
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
  }
};

double derivest(double(*fun)(double), double x0, DerivestParams par, double* err, double* finaldelta) {
  const double step_ratio = 2.0000001;  // DERIVESTcpp is hardcoded to this specific step_ratio in several places.
  const double h = std::max(x0, 0.02);  // same as nominal_step

  static const int ndel = 26;
  double delta[ndel];
  for (int i = 0; i < ndel; i++) delta[i] = par.max_step*pow(step_ratio, -i);
  
  // generate finite differencing rule in advance.
  // The rule is for a nominal unit step size, and will
  // be scaled later to reflect the local step size.
  double fdarule[10]; int nfda = 0;
  get_fdarule(par.derivative_order, par.method_order, par.style, fdarule, &nfda);
  
  // will we need fun(x0)?
  const double f_x0 = (par.derivative_order % 2 == 0 || par.style != DerivestStyle_Central) ? fun(x0) : 0.0;

  double der = 0.0, errest = 0.0;
  if (finaldelta != nullptr) *finaldelta = 0.0;

  double f_del[ndel];
  if (par.style == DerivestStyle_Central) {
    // A central rule, so we will need to evaluate symmetrically around x0i.
    double f_plusdel[ndel], f_minusdel[ndel];
    for (int i = 0; i < ndel; i++) f_plusdel[i] = fun(x0 + h*delta[i]);
    for (int i = 0; i < ndel; i++) f_minusdel[i] = fun(x0 - h*delta[i]);
    if (par.derivative_order == 1 || par.derivative_order == 3) { for (int i = 0; i < ndel; i++) f_del[i] = (f_plusdel[i] - f_minusdel[i]) / 2.0; } // odd transformation
    else { for (int i = 0; i < ndel; i++) f_del[i] = (f_plusdel[i] + f_minusdel[i]) / 2.0 - f_x0; }
  } else if (par.style == DerivestStyle_Forward) {
    for (int i = 0; i < ndel; i++) f_del[i] = fun(x0 + h*delta[i] - f_x0);  // forward rule; drop off the constant only
  } else if (par.style == DerivestStyle_Backward) {
    for (int i = 0; i < ndel; i++) f_del[i] = fun(x0 - h*delta[i] - f_x0);  // backward rule; drop off the constant only
  }

  // Apply the finite difference rule at each delta, scalingas appropriate for delta and the requested DerivativeOrder.
  // First, decide how many of these estimates we will end up with.
  const int ne = ndel + 1 - nfda - par.romberg_terms;

  // Form the initial derivative estimates from the chosen finite difference method.
  // der_init = vec2mat(f_del, ne, nfda)*fdarule.'
  double der_init[ndel];  
  for (int i = 0; i < ne; i++) {
    der_init[i] = 0;
    for (int j = 0; j < nfda; j++) der_init[i] += f_del[i + j] * fdarule[j];
  }

  // scale to reflect the local delta (" der_init = der_init(:). / (h*delta(1:ne)).^par.DerivativeOrder " )
  for (int i = 0; i < ne; i++) der_init[i] = der_init[i] / pow(h*delta[i], par.derivative_order);

  // Each approximation that results is an approximation of order par.DerivativeOrder to the desired derivative.
  // Additional(higher order, even or odd) terms in the Taylor series also remain. Use a generalized (multi-term)
  // Romberg extrapolation to improve these estimates.
  int rombexpon[3] = { -1,-1,-1 };
  if (par.style == DerivestStyle_Central) for (int i = 0; i < par.romberg_terms; i++) rombexpon[i] = 2 * i + par.method_order - 2;
  else for (int i = 0; i < par.romberg_terms; i++) rombexpon[i] = i + par.method_order - 1;

  
  // for (int i = 0; i < ne; i++) printf("%.6f\t", der_init[i]); printf("\n");
  return 0.0;
}

int main() {
  derivest(exp, 0, DerivestParams(), nullptr, nullptr);
  double* ar = get_fdarule_test(0);
  for (int i = 0; i < 3; i++) printf("%.3f ", ar[i]);

  ar = get_fdarule_test(1);
  for (int i = 0; i < 3; i++) printf("%.3f ", ar[i]);
  return 0;
}