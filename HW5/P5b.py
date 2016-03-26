#####################################################################################################
# Milos Atz
# NE155 Homework 5
#####################################################################################################
import math
import numpy as np
import scipy
#####################################################################################################
# Problem 5b
# First, build the system using the method performed in Problem 1
def matrix_build(n):
	a=[-1]*int(n-1)
	b=[4]*int(n)
	c=[-1]*int(n-1)
	A=np.matrix(np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1))
	return(A)
def b_build(n):
	b=np.zeros(n)
	for i in range(0, n):
		b[i]=100
	b=np.transpose(np.matrix(b))
	return(b)
#####################################################################################################
# Then, set up the vector of omegas to test.
w=[None]*101
w_min=0
w_max=2
w_step=w_max/float(len(w)-1)
for i in range(0,len(w)):
	w[i]=w_min+i*w_step
#####################################################################################################
# Initialize the D, L, and U matrices from which to develop the P matrix.
A=matrix_build(5)
b=b_build(5)
D=np.diag(np.diag(A))
L=np.diag(np.diag(A,-1),-1)
U=np.diag(np.diag(A,1),1)
#####################################################################################################
# To find the optimum omega, we find the spectral radius for each P matrix. The one that yields the P matrix with the lowest spectral radius is that which is the optimum.
specRad=[None]*101
for i in range(0, len(w)):
	P_sor=np.linalg.inv(D+w[i]*L)*(D*(1-w[i])-w[i]*U)
	specRad[i]=max(abs(np.linalg.eigvals(P_sor)))
w_opt=w[specRad.index(min(specRad))]
print('optimum omega = ' + str(w_opt))
#####################################################################################################
import matplotlib.pyplot as plt
fig=plt.plot(w, specRad, label='spectral radius')
plt.ylabel('spectral radius')
plt.xlabel('omega')
plt.savefig('p5b.png', bbox_inches='tight')
plt.clf()