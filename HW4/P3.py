#################################################################################
# Milos Atz
# NE155 Homework 3
#################################################################################
import math
import numpy
import scipy
##########################################################################
# Problem 3: Solve the following systems of linear equations with:
b=numpy.matrix('20; -7; 4; 6')

# System 1
sys1=numpy.matrix('4 -1 2 3; 0 -2 7 -4; 0 0 6 5; 0 0 0 3')
x1 = numpy.linalg.solve(sys1, b)

# System 2
#sys2=numpy.matrix('4 -1 2 3; 0 0 7 -4; 0 0 6 5; 0 0 0 3')
#x2 = numpy.linalg.solve(sys2, b)
# Throws error

# System 3
#sys3=numpy.matrix('4 -1 2 3; 0 0 7 0; 0 0 6 5; 0 0 0 3')
#x3 = numpy.linalg.solve(sys3, b)
# Throws error
