#! /usr/bin/env python3
import scipy.integrate
import math

url = '''https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html'''

print("See full documentation: ", url)


def integrand(x, omega):
    return math.sin(omega*x)


a = 0.0
b = 2*math.pi

omega = 1.7


# Because our function will end up having very small values, it may be useful to set the 'absolute error' goal to zero, and use the relative error

integral, err = scipy.integrate.quad(
    integrand, a, b, args=(omega,), epsabs=0.0, epsrel=1.0e-3)

print(integral)
