#! /usr/bin/env python2
import scipy.optimize
import numpy as np
import math

url = '''https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html'''

print("See full documentation: ", url)


def f(x):
    return math.sin(x)

# Want to solve f(x)=f0, for x
# ==> f(x)-f0=0
# Need to give a range, where we know solution lies inside range


f0 = 0.5


def function(x):
    return f(x) - f0


print(function(0), function(math.pi/2))

x_value = scipy.optimize.brentq(function, 0.0, math.pi/2)

print("x value: ", x_value, ", function: f(x)=", f(x_value))
