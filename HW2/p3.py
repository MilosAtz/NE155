#################################################################################
# Milos Atz
# NE155 Homework 2
#################################################################################
import math
import numpy
import scipy
import types
import matplotlib.pyplot as plt
#################################################################################
# Problem 3
# We have the following data: x = [1, 2, 3, 4, 5, 6], f(x) = [1, 2, 15, 12, 7, 3].
# a) Using built-in Python functions, interpolate this data using (i) piecewise linear interpolation, (ii) lagrange polynomial interpolation, and (iii) spline interpolation
#################################################################################
xStart=0.75
xEnd=6.25
xIncr=0.05
nPoints=int((xEnd-xStart)/xIncr)
xPlot=[None]*(nPoints+1)
for i in range(0,len(xPlot)):
	xPlot[i]=xStart+i*xIncr

x = [1, 2, 3, 4, 5, 6]
fx = [1, 2, 15, 12, 7, 3]
#################################################################################
# (i) Piecewise linear interpolation:
# numpy.interp performs one-dimensional linear interpolation; returns the one-dimensional piecewise linear interpolant to a function with given values at discrete data-points.

y_i=numpy.interp(xPlot, x, fx)

figData=plt.plot(x, fx, 'ro', label='Data points')
figInterp=plt.plot(xPlot, y_i, label='Piecewise linear interpolated values')
plt.ylabel('f(x)')
plt.xlabel('x')
plt.legend(loc=4)
plt.savefig('p3a_i.png', bbox_inches='tight')
plt.clf()

#################################################################################
# (ii) Lagrange polynomial interpolation
# scipy.interpolate.lagrange(x, w)returns a Lagrange interpolating polynomial, given two 1-D arrays x and w, returns the Lagrange interpolating polynomial through the points (x, w). Warning: This implementation is numerically unstable. Do not expect to be able to use more than about 20 points even if they are chosen optimally.

from scipy.interpolate import lagrange

y_ii=lagrange(x,fx)
yVals_ii=[None]*len(xPlot)
for i in range(0,len(xPlot)):
	yVals_ii[i]=y_ii(xPlot[i])

figData=plt.plot(x, fx, 'ro', label='Data points')
figInterp=plt.plot(xPlot, yVals_ii, label='Lagrange polynomial interpolated values')
plt.ylabel('f(x)')
plt.xlabel('x')
plt.legend(loc=4)
plt.savefig('p3a_ii.png', bbox_inches='tight')
plt.clf()

#################################################################################
# (iii) Spline interpolation (piecewise polynomial interpolation)
from scipy.interpolate import splev, splrep
y_iii = splrep(x, fx)
yVals_iii = splev(xPlot, y_iii)
figData=plt.plot(x, fx, 'ro', label='Data points')
figInterp=plt.plot(xPlot, yVals_iii, label='Spline interpolated values')
plt.ylabel('f(x)')
plt.xlabel('x')
plt.legend(loc=4)
plt.savefig('p3a_iii.png', bbox_inches='tight')
plt.clf()

#################################################################################
# b) Briefly discuss the differences between the interpolations

# The linear interpolation is the simplest of the three. It produces straight lines between the data points. Thus, the resulting interpolation is not smooth and assumes that the behavior of the function that we're interpolating is between the points. We assume that the maximum and minimum are defined by the data points provided. Beyond the first and last points, the interpolated function is flat.

# The lagrange interpolation adds smoothness to the interpolated result, but there are some interesting changes in inflection. While the linear interpolation produced what looked like a rough sketch of a normal bell curve (2 inflection points), the lagrange interpolation has 4 inflection points. Due to the way the points are arranged, the interpolation predicts a large dip between the first two points and a change in inflection between the last two points.

# The spline interpolation is an even smoother version of the lagrange interpolation. The dip between the first two points is not as deep and the change in inflection between the last two points is gone, providing a smoother, more toned result. The general tendency of the lagrange interpolation is maintained, and both the spline and lagrange interpolations are very different from the linear interpolation.





