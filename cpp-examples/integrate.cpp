#include <cassert>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <vector>

/*
GSL provides many functions to perform integrals.
https://www.gnu.org/software/gsl/doc/html/integration.html

Many use some form of adaptive Gaussian quadrature, which are typically the
most accurate methods to use.
Each is best for different cases

All the functions take as input a pointer to a 'gsl_function'
Which is GSLs way to evaluate a function of any number of parameters

Some common examples:

The "workhorse" is the qag() and qags() functions
gsl_integration_qag()
 * Integrates from [a, b]

gsl_integration_qags()
 * Integrates from [a, b]
 * Specialised for integrands with singularities

Special functions also given for (semi-) infinite intervals

gsl_integration_qagi()
 * Integrates from [-infty, +infty]

gsl_integration_qagiu()
 * Integrates from [a, +infty]

gsl_integration_qagil()
 * Integrates from [-infty, b]

*/

//! this is the function we will integrate. We integrate over x, but it also
//! depends on a -- this is to demonstrate how to deal with extra parameters
double function(double x, double y, double z) {
  return std::exp(-y * x) * std::cos(x / z);
}

//! gsl_function must be of this form
double function_gsl_form(double x, void *params) {
  // "unpack" the void* into the actual parameters
  // (this can be done in fewer steps, but I want to make it clearer)
  assert(params != nullptr);
  const auto array = static_cast<double *>(params);
  const auto a = array[0];
  const auto b = array[1];
  return function(x, a, b);
}

double indef_integral(double x, double y, double z) {
  return std::exp(-y * x) * (-y * z * std::cos(x / z) + std::sin(x / z)) /
         (1.0 + y * y * z * z);
}

int main() {

  // Function parameters: simple case: y=z=1
  double y = 1.0;
  double z = 1.0;
  // Integration range: [-2*pi, 2*pi]
  double a = -2.0 * M_PI;
  double b = 2.0 * M_PI;

  // set the relativa and absolute error goals:
  double abs_err = 1.0e-4;
  double rel_err = 1.0e-4;
  // nb: in assignment, our function will end up having very small values, it
  // may be useful to set the 'absolute error' goal to zero, and use the
  // relative error

  // set up the inputs (gsl_function) GSL integrator
  std::vector<double> params{y, z};
  gsl_function f_gsl;
  f_gsl.function = function_gsl_form;
  f_gsl.params = params.data();

  // set up the options for GSL integrator, and allocate workspace memory
  const unsigned long max_num_subintvls = 1000;
  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);

  // exact result integrated from [-2*pi, 2*pi] is Sinh[2*pi] =~267.745..
  double result{0.0};
  gsl_integration_qag(&f_gsl, a, b, abs_err, rel_err, max_num_subintvls,
                      GSL_INTEG_GAUSS15, gsl_int_wrk, &result, &abs_err);

  // free workspace memory
  gsl_integration_workspace_free(gsl_int_wrk);

  // expected result:
  const auto expected = indef_integral(b, y, z) - indef_integral(a, y, z);
  const auto error = result - expected;

  std::cout << result << " " << expected << " " << error << "\n";
  // Notice the result is *much* more accurate than the accuracy we requested.
  // This is simply because this is an easy function to integrate. We won't be
  // so lucky in general
}