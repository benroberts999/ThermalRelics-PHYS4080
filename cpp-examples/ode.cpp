#include <cassert>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <vector>

/*!
Requires gsl libraries ('apt install libgsl-dev' or 'brew install gsl')
Compile, e.g., with: g++ ode.cpp -lgsl -lblas

This is an c++ example program that solves the ODE:
  y'(t) = y(t) + a*Sin(t)
with initial condition:
  y(0) = y0


The analytic solution is:
y(t) = (a/2) * [exp(t){1 + 2(y0/a)} - cos(t) - sin(t)]

It uses the GSL (GNU Scientific Library).
see: https://www.gnu.org/software/gsl/doc/html/ode-initval.html

*/

//------------------------------------------------------------------------------

//! The derivative we wish to solve the ODE for, dy/dt. First argument is t,
//! second is y itself followed by any other parameters (in this case, just a)
double dydt_ode(double t, double y, double a) {
  // y'(t) = y(t) + a*Sin(t)
  return y + a * std::sin(t);
}

//! Now we must define a "wrapper" function, using the language expected by the
//! GSL library. This will be passed into GSL library
int ode_function_gsl(double t, const double y[], double dydt[], void *params) {

  // define ODE function in language of GSL:
  // It must have the following signature:
  // int (* function) (double t, const double y[], double dydt[], void * params)

  // y[] and dydt[] are _arrays_, since general ODE may be of more than 1
  // variable. In our case, these are simple 1D arrays since we have a 1D
  // problem.

  // The void pointer 'params' allows us to pass arbitrary parameters to the
  // dydt function. Usually, this is a pointer to an array of doubles, but can
  // in general be anything (e.g., a Struct)
  // This function must fill dydt[i] using y[i] and params (for ith equation in
  // system) - for us, 1D system
  // see: https://www.gnu.org/software/gsl/doc/html/ode-initval.html

  // // Interpret params (GSL needs as a void pointer)
  // const double *const params_array = static_cast<double *>(params);
  // // Extract parameters (for us, single parameter, a)
  // const double a = params_array[0];
  const double *a_ptr = static_cast<double *>(params);
  assert((a_ptr != nullptr) && "in ode_function_gsl, params may not be null");
  // fill the dydt array
  dydt[0] = dydt_ode(t, y[0], *a_ptr);
  return GSL_SUCCESS;
}

//------------------------------------------------------------------------------
// We define a "wrapper" function, which calls GSL and solves the ode
// We use templates to allow passing of any type of function, and any set of
// paramers
// Note: 'Function' must have signature expected by GSL
// see: https://www.gnu.org/software/gsl/doc/html/ode-initval.html
template <typename Function, typename Parameters>
std::vector<double> solve_ode(const std::vector<double> &t_list,
                              double y_initial, Function func,
                              Parameters params) {

  // define array to store solution:
  std::vector<double> y_solution;
  y_solution.reserve(t_list.size());

  // Set up the ODE system for GSL
  static const int ode_dim = 1;  // single first-order ODE
  const auto Jacobian = nullptr; // 1D ODE, so do not require Jacobian
  gsl_odeiv2_system sys = {func, Jacobian, ode_dim, &params};

  // Some settings for the solver:
  const auto ode_method = gsl_odeiv2_step_rkf45; // method (RK4)
  const double hstart = 1.0e-6;                  // initial step size
  const double epsabs = 1.0e-3;                  // absolute error goal
  const double epsrel = 1.0e-3;                  // relative error goal
  gsl_odeiv2_driver *ode_solver =
      gsl_odeiv2_driver_alloc_y_new(&sys, ode_method, hstart, epsabs, epsrel);

  // Set the [t,y(t)] with initial values
  double y[ode_dim] = {y_initial};
  double t = t_list.front();

  // loop through each t for which we want the solution:
  for (auto next_t : t_list) {

    // Drive the ODE system from t -> t_i = t + dt
    // (i.e., find value of y(t) at new t=t+dt, using value at previous t)
    const auto status = gsl_odeiv2_driver_apply(ode_solver, &t, next_t, y);
    if (status != GSL_SUCCESS) {
      std::cout << "Warning: didn't converge at t=" << t << ", y=" << y[0]
                << "\n";
    }
    // store solution:
    y_solution.push_back(y[0]);
  }

  // Free memory ascosciated with driver:
  gsl_odeiv2_driver_free(ode_solver);

  return y_solution;
}

//------------------------------------------------------------------------------
int main() {

  // a is parameter of ODE;
  const double a = 3.0;
  // Initial condition:
  const double y0 = 1.0;

  // set up "grid" of t values to solve y(t)
  const double t0 = 0.0;
  const double tmax = 4.0;
  const int n_steps = 100;
  const double dt = (tmax - t0) / (n_steps + 1);
  // store list of t's to solve for:
  std::vector<double> t_list;
  for (int i = 0; i < n_steps; ++i) {
    t_list.push_back(t0 + i * dt);
  }

  // solve the ode:
  const auto y_sol = solve_ode(t_list, y0, ode_function_gsl, a);

  // write solution (and error) to screen
  for (std::size_t i = 0; i < t_list.size(); ++i) {
    const auto t = t_list.at(i);
    const auto y_numerical = y_sol.at(i);
    const auto y_exact = (a / 2.0) * (std::exp(t) * (1.0 + 2.0 * (y0 / a)) -
                                      std::cos(t) - std::sin(t));
    std::cout << t << " " << y_numerical << " error:" << y_numerical - y_exact
              << "\n";
  }
}
