#! /usr/bin/env python3
"""
Example, which solves the differential equation:
dy/dt = a * t * cos(b*t) * y(t)
This has simple analytic solution:
y = c * Exp[a(cos(b*t)/b^2 + t * sin(b*t)/b)]
We'll solve with y[0]=1, so c = Exp[-a]
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint

url2 = """https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html"""

print("See full documentation: ", url2)


def dy_dt(t, y, a1, a2):
    """Must have form (t, y, params....)"""
    return a1 * t * np.cos(a2 * t) * y


# Here we use linear 'grid', but often logarithmic is better!
t_list = np.linspace(0.0, 15.0, 1000)


def solve_dYdx(t_list, y0, a1, a2):
    """We don't need to wrap this in a function, but can be convenient...
    Particularly if there are other things we need to calculate first, for example y[t0]
    """
    return odeint(dy_dt, y0, t_list, args=(a1, a2), tfirst=True)


y0 = 1.0
a = 1.5
b = 6.3
sol = solve_dYdx(t_list, y0, a, b)


# Analytic solution
def y_analytic(t, y0, a, b):
    return (
        y0
        * np.exp(-a / b**2)
        * np.exp(a * (t * np.sin(b * t) / b + np.cos(b * t) / b**2))
    )


exact_sol = y_analytic(t_list, y0, a, b)

plt.plot(t_list, sol, ".-", label="Numerical solution (odeint)")
plt.plot(t_list, exact_sol, label="Exact solution")
plt.legend(loc="best")
plt.ylabel("y")
plt.xlabel("t")
plt.grid()
plt.show()
