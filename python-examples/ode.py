#! /usr/bin/env python3
'''
Example, which solves the simple differential equation:
dy/dt =  y'(t) = y(t) + a*Sin(t)
with initial condition:
  y(0) = y0
This has analytic solution:
y(t) = (a/2) * [exp(t){1 + 2(y0/a)} - cos(t) - sin(t)]
'''

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

url1 = '''https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html'''
url2 = '''https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html'''

print("See full documentation: ", url1, "\n", url2)


def dydt(t, y, a):
    return y + a * math.sin(t)


t0 = 0.0
tmax = 4.0

#This is y(t0)
y0 = 1.0
a = 3.0

print("Values used: y0={}, a={}".format(y0, a))

# We pass the "other" arguments for the function (dydt) in as a list
# 't' should be the first argument of the dy/dt(t,...) function

t_list = np.linspace(t0, tmax, 100)
sol = solve_ivp(dydt, [t0, tmax], [y0], t_eval=t_list,
                args=[a], rtol=1.0e-3, atol=0.0)


sol2 = odeint(dydt, y0, t_list, args=(a,), tfirst=True)


def exact_sol(t):
    return (a / 2.0) * (np.exp(t) * (1.0 + 2.0 * (y0 / a)) -
                        np.cos(t) - np.sin(t))


plt.plot(sol.t, sol.y[0], label='Numerical solution (solve_ivp)')
plt.plot(t_list, sol2, label='Numerical solution (odeint)')
plt.plot(sol.t, exact_sol(sol.t), label='Exact solution')
plt.legend(loc='best')
plt.ylabel('y')
plt.xlabel('t')
plt.grid()
plt.show()
