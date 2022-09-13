#! /usr/bin/env python3
import scipy.integrate
import math

url = '''https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html'''

print("See full documentation: ", url)


# We wish the integrate this function of x. Function also contains other parameters
def integrand(x, y, z):
    return math.exp(-y * x) * math.cos(x / z)


# This is the exact (indefinite) integral
def indef_integral(x,  y,  z):
    return math.exp(-y * x) * z * (-y * z * math.cos(x / z)
                                   + math.sin(x / z)) / (1.0 + y * y * z * z)


# We integrate [-2pi, 2pi]
a = -2*math.pi
b = 2*math.pi

# simple case: other parameters y=z=1
y = 1.0
z = 1.0


# Because our function will end up having very small values, it may be useful to set the 'absolute error' goal to zero, and use the relative error

integral, error_estimate = scipy.integrate.quad(
    integrand, a, b, args=(y, z), epsabs=1.0e-6, epsrel=1.0e-6)

print(integral, error_estimate)

expected = indef_integral(b, y, z) - indef_integral(a, y, z)
error = integral - expected

print(expected, error)
