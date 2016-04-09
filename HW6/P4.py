#################################################################################
# Milos Atz
# NE155 Homework 6
#################################################################################
import math
import numpy as np
import scipy.linalg
import math
#####################################################################################################
# Problem 4
# Numerically solve the eigenvalue form of the diffusion equation with the same BCs as in Problem 2. Use FDM for the discretization of the spatial variable; use Power Iteration to find the dominant eigenvalue and corresponding eigenvector. Note that you may need to normalize the solution vector (should at least normalize the initial guess). Use SOR or GS method to complete the solve portion of the algorithm. Use absolute error tolerance to check for convergence.
#####################################################################################################
# Define the GS solver used in the problem script; this is the same as what was used in HW5.
def gs_solver(A, b, tol=1e-6):
	if(min(np.linalg.eigvals(A)<0)):
		sys.exit('A is not positive definite')
	if((A.transpose() != A).all()):
		sys.exit('A is not symmetric')
	if(b.size!=A.shape[0] or b.size!=A.shape[1]):
		sys.exit('dimensions of A and b do not agree')
	n=b.size
	x_old=np.transpose(np.matrix(np.zeros(n)))
	D=np.diag(np.diag(A))
	L=np.diag(np.diag(A,-1),-1)
	U=np.diag(np.diag(A,1),1)
	DL_inv=np.linalg.inv(D+L)
	conv=1
	counter=0
	while(conv>tol):
		x_new=DL_inv*(-U*x_old+b)
		# print(x_new)
		conv=np.linalg.norm(x_new-x_old)
		x_old=x_new
		counter=counter+1
	# print('counter= '+str(counter))
	# print('absolute error = '+str(conv))
	return(x_new)
#####################################################################################################
a = 4.0			# cm
D = 1.0			# cm
sig_a = 0.7		# 1/cm
vsig_f = 0.6	#1/cm
h = 0.1			# cm
#####################################################################################################
# First, we determine the number of cells and points.
n_cell = int((a-(-1*a))/h)
n_points = n_cell+1
#####################################################################################################
# We need initial values for k(0) and phi(0) for i=0,...,n-1; we normalize phi0 = phi0/norm(phi0). Let's assume that k(0)=1 and phi(0)=1 for i=0,...,n-1.
k = 0.98
phi = np.transpose(np.matrix(np.ones(n_cell-1)))
phi = phi/np.linalg.norm(phi)
#####################################################################################################
# We compute A in the same way as done for Problem 2. A is made up out of the coefficients for flux; A is a tridiagonal matrix. The inputs a, b, and c allow for the input of those coefficients.
A_a = [-1]*int(n_cell-2)
A_b = [2+(sig_a*h**2/D)]*int(n_cell-1)
A_c = [-1]*int(n_cell-2)
A = np.matrix(np.diag(A_a, -1)+np.diag(A_b, 0)+np.diag(A_c, 1))
#####################################################################################################
# The initial fission source is calculated as nu*sig_F_i,i*phi(0)_i. nu*sig_F is constant, so we just multiply that by the initial guess for phi.
Q = (h**2*vsig_f/D)*phi
#####################################################################################################
# In order to iterate until convergence, we set up a convergence criterion, defined by the error of the solution relative to our desired tolerance. We have two criteria to set up - one for the eigenvector (phi) and one for the eigenvalue (k). The absolute error will be calculated as k(n)-k(n-1) and ||Ax-b|| for phi.
error_k = 1
error_phi = 1
tol = 1e-4
counter = 0
#####################################################################################################
# Within the while loop, we will iterate to solve A*phi = (1/k)*Q using the GS method, which has it's own iteration loop.
while(error_k > tol or error_phi > tol):
	counter = counter+1
	b = (1/k)*Q
	# Use the GS solver to solve for phi(m)
	phi_new = gs_solver(A,b)
	# Compute the next fission source
	Q_new = (h**2*vsig_f/D)*phi_new
	# Compute the next eigenvalue
	k_new = k*sum(Q_new)/sum(Q)
	# Check for convergence
	error_k = abs(k_new-k)
	error_phi = np.linalg.norm(phi_new-phi)
	k=float(k_new)
	# print(k)
	phi=phi_new
	Q=Q_new
print('k = '+str(k))
print('number of power iterations = '+str(counter))
print('error in phi = '+str(float(error_phi)))
print('error in k = '+str(float(error_k)))
phi=np.transpose(np.insert(phi, 0, 0.0))
phi=np.transpose(np.insert(phi, n_points-1, 0.0))
#####################################################################################################
# Plot the eigenvector, phi, from x = -a to a, Report the eigenvalue, k, and the number of power iterations required for convergence.
import matplotlib.pyplot as plt
x_vals=[0]*n_points
for i in range(0,n_points):
	x_vals[i]=-1*a+h*i
x_vals=np.transpose(np.matrix(x_vals))
fig=plt.plot(x_vals, phi, 'rx', label='finite difference')
plt.xlim([-a-0.5,a+0.5])
plt.ylabel('phi(x)')
plt.xlabel('x')
plt.figtext(0.14, 0.86, 'k = '+str(round(k,5)))
plt.figtext(0.14, 0.82, 'power iterations = '+str(counter))
plt.figtext(0.14, 0.78, 'error in phi = '+str(round(float(error_phi),9)))
plt.figtext(0.14, 0.74, 'error in k = '+str(round(float(error_k),9)))
#plt.show()
plt.savefig('p4_flux.png', bbox_inches='tight')
plt.clf()


