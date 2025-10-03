#! /usr/bin/env python3
import scipy.optimize
import numpy as np

url = """https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html"""

print("See full documentation: ", url)


# Function we want to solve: sin(x) = f0
def f1(x, a, b):
    return a * np.sin(b * x)


# Want to solve f(x)=f0, for x
# ==> f(x)-f0=0
# Need to give a range, where we know solution lies inside range
f0 = 0.5


# Define the function whose root we seek: f(x) - f0 = 0
def f1_root(x, a, b):
    return f1(x, a, b) - f0


# Initial guess (instead of bracket, required by fsolve)
x_guess = 0.5

a, b = 1.0, 3.0

# Solve
x_value = scipy.optimize.fsolve(f1_root, x_guess, args=(a, b))[0]

print(f"x value: {x_value:.6f}, f(x) = {f1(x_value,a,b):.6f}")
