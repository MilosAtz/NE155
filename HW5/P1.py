#################################################################################
# Milos Atz
# NE155 Homework 5
#################################################################################
import math
import numpy as np
import scipy.linalg
##########################################################################
# Problem 1
# a) Use built in Python commands to construct A and b.
def matrix_build(n):
	a=[-1]*int(n-1)
	b=[2]*int(n)
	c=[-1]*int(n-1)
	A=np.matrix(np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1))
	return(A)
def b_build(n):
	b=np.zeros(n)
	for i in range(0, n):
		b[i]=i
	b=np.transpose(np.matrix(b))
	return(b)
##########################################################################
# b) What is the condition number of A?
n=100
A=matrix_build(n)
conditionN=np.linalg.cond(A)
##########################################################################
# c) Solve this problem by explicitly inverting A and multiplying b.
b=b_build(n)
A_inv=np.linalg.inv(A)
soln_c=A_inv*b
##########################################################################
# d) Solve this problem by explicitly inverting A and multiplying b.
soln_d=scipy.linalg.solve(A, b)
##########################################################################
# plot the solutions for (c) and (d)
import matplotlib.pyplot as plt
fig_c=plt.plot(b, soln_c, label='solution (c)')
fig_d=plt.plot(b, soln_d, 'rx', label='solution (d)')
plt.ylabel('solution (x)')
plt.xlabel('b')
plt.legend(loc=2)
#plt.show()
plt.savefig('p1.png', bbox_inches='tight')
plt.clf()
