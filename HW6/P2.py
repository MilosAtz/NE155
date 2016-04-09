#################################################################################
# Milos Atz
# NE155 Homework 6
#################################################################################
import math
import numpy as np
import scipy.linalg
import math
#####################################################################################################
# Problem 2
# Numerical solve the fixed-source diffusion equation as described in Question 1 using the finite difference method for discretization of the spatial variable and Gaussian elimination (aka the Thomas algorithm; note that you will have a tridiagonal system to solve) for solving the system of linear algebraic equations.
#####################################################################################################
a = 4.0		# cm
D = 1.0		# cm
sig_a = 0.2	# 1/cm
S = 8.0		# n/(cm^3*s)
h = 0.1		# cm
#####################################################################################################
# First, we determine the number of cells and points.
n_cell = int((a-(-1*a))/h)
n_points = n_cell+1
#####################################################################################################
# We want to set up the system Ax=b, where A is a tridiagonal matrix and x contains the flux at each point. Because S is constant, the vector b will be constant as well.
b=np.zeros(n_cell-1)
for i in range(0, n_cell-1):
	b[i]=S*h**2/D
b=np.transpose(np.matrix(b))
#####################################################################################################
# A is made up out of the coefficients for flux; A is a tridiagonal matrix. The inputs a, b, and c allow for the input of those coefficients.
A_a=[-1]*int(n_cell-2)
A_b=[2+(sig_a*h**2/D)]*int(n_cell-1)
A_c=[-1]*int(n_cell-2)
A=np.matrix(np.diag(A_a, -1) + np.diag(A_b, 0) + np.diag(A_c, 1))
#####################################################################################################
# Utilize the Thomas method to solve the system of equations
phi = [0]*n_points
for i in range(1,n_cell-1):
	A[i,i] = A[i,i]-(A[i,i-1]/A[i-1,i-1])*A[i-1,i]
	b[i]=b[i]-(A[i,i-1]/A[i-1,i-1])*b[i-1]

phi[n_cell-1]=b[n_cell-2]/A[n_cell-2,n_cell-2]
for i in range(n_cell-3, -1, -1):
	phi[i+1]=(b[i]-A[i,i+1]*phi[i+2])/A[i,i]
#####################################################################################################
# Plot the solution from x = -a to x = a. Compare your answer to one of your solutions from Q1.
import matplotlib.pyplot as plt
x_vals=[0]*n_points
for i in range(0,n_points):
	x_vals[i]=-1*a+h*i
# The analytical solution from 1b should be plotted to compare the numerical result
phi_a = [0]*n_points
L=math.sqrt(D/sig_a)
for i in range(0,n_points):
	phi_a[i]=(-S/sig_a)*((math.exp(-x_vals[i]/L)+math.exp(x_vals[i]/L))/(math.exp(-a/L)+math.exp(a/L))-1)
fig=plt.plot(x_vals, phi, 'rx', label='finite difference')
fig=plt.plot(x_vals, phi_a, label='analytical solution')
plt.xlim([-a-0.5,a+0.5])
plt.ylim([-0.5, int(math.ceil(max(phi_a)/10.0))*10+5])
plt.ylabel('phi(x)')
plt.xlabel('x')
plt.legend(loc=2)
#plt.show()
plt.savefig('p2_flux.png', bbox_inches='tight')
plt.clf()
# Finally, compare the results by plotting the difference between them.
diff=[0]*n_points
for i in range(0,n_points):
	diff[i]=phi_a[i]-phi[i]
fig=plt.plot(x_vals, diff, label='difference between analytical and numerical solutions')
plt.ylabel('difference')
plt.xlabel('x')
#plt.show()
plt.savefig('p2_difference.png', bbox_inches='tight')
plt.clf()
#####################################################################################################
# The result of the numerical solution is very close to that of the analytical solution plotted over the same space. The largest discrepancy occurs at x=0, where the flux is the highest.
#####################################################################################################

