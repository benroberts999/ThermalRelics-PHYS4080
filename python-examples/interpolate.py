#! /usr/bin/env python3
import scipy.interpolate
import numpy as np

url = '''https://docs.scipy.org/doc/scipy/reference/interpolate.html'''

print("See full documentation: ", url)

x_data = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])

y_data = np.array([0.0, 1.0, 4.0, 9.0, 16.0, 25.0])

# interpolator
y_interp = scipy.interpolate.CubicSpline(
    x_data, y_data, extrapolate=True)

# Note: exptrapolating using interpolation functions is often a VERY bad idea - always check if this is reasonable for your problem

x_list = [1.0, 1.5, 4.5, 5.0, 6.0]
for x in x_list:
    print("x=", x, " y_interp(x)=", y_interp(x))
