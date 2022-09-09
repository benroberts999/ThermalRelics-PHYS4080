#! /usr/bin/env python3
'''
Example, which solves the "trivial" differential equation:
dy/dt = a * cos(t)
This has simple analytic solution:
y = y0 + a * sin(t)
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
    return a*math.cos(t)


t0 = 0.0
tmax = 2*math.pi

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


exact_sol = 1.0 + a*np.sin(sol.t)

plt.plot(sol.t, sol.y[0], label='Numerical solution (solve_ivp)')
plt.plot(t_list, sol2, label='Numerical solution (odeint)')
plt.plot(sol.t, exact_sol, label='Exact solution')
plt.legend(loc='best')
plt.ylabel('y')
plt.xlabel('t')
plt.grid()
plt.show()
