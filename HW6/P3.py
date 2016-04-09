#################################################################################
# Milos Atz
# NE155 Homework 6
#################################################################################
import math
import numpy as np
import scipy.linalg
import math
import matplotlib.pyplot as plt
#####################################################################################################
# Problem 3
# Investigate how well your numerical solution approximates the analytical solution by compute phi_i for various constant mesh sizes. For each mesh length calculate the relative error between your numerical and analytical solutions. Plot the maximum relative error as a function of total number of meshes for each case. What can you conclude about the relationship between maximum error and total number of meshes?
#####################################################################################################
# Turn everything from Problem 2 in to a function!
def flux_fdm(h,a,D,sig_a,S):
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
	return(phi)
#####################################################################################################
# We now have a function that will return phi for variable parameter values. We want to find the maximum error between our FDM result for flux and the analytical solution. Plot the solution from x = -a to x = a. Compare your answer to one of your solutions from Q1.
a = 4.0		# cm
D = 1.0		# cm
sig_a = 0.2	# 1/cm
S = 8.0		# n/(cm^3*s)
L=math.sqrt(D/sig_a)
h = [1.0, 0.5, 0.1, 0.05, 0.01]		# cm
n_points = [0]*len(h)
#####################################################################################################
# Evaluate the FDM result for various values of h. Find the maximum relative error between the analytical and FDM result and plot that against h.
RE_max=[0]*len(h)
for j in range(0, len(h)):
	n_points[j] = int((a-(-1*a))/h[j])+1
	# Evaluate the FDM result of phi for the given spacing h
	phi=flux_fdm(h[j],a,D,sig_a,S)
# Create a vector of x values on which to solve the exact analytical equation.
	x_vals=[0]*n_points[j]
	for i in range(0,n_points[j]):
		x_vals[i]=-1*a+h[j]*i
# The analytical solution from 1b should be plotted to compare the numerical result
	phi_a = [0]*n_points[j]
	RE = [0]*n_points[j]
	for i in range(1,n_points[j]-1):
		phi_a[i]=(-S/sig_a)*((math.exp(-x_vals[i]/L)+math.exp(x_vals[i]/L))/(math.exp(-a/L)+math.exp(a/L))-1)
		RE=abs(phi_a[i]-phi[i])/phi_a[i]
	RE_max[j]=max(RE)
#####################################################################################################
# Finally, plot the maximum relative error as a function of h
fig=plt.plot(n_points, np.squeeze(RE_max))
plt.ylabel('Maximum relative error')
plt.xlabel('Number of Points')
plt.yscale('log')
#plt.show()
plt.savefig('p3_re.png', bbox_inches='tight')
plt.clf()
#####################################################################################################
# 'The maximum error decreases exponentially as the number of meshes is increased! More meshes means drastically better results!'
print('What can you conclude about the relationship between the maximum error and the total number of meshes?'+'\n'+'The maximum error decreases exponentially as the number of meshes is increased! More meshes means drastically better results!')
#####################################################################################################





