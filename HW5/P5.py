#####################################################################################################
# Milos Atz
# NE155 Homework 5
#####################################################################################################
import math
import numpy as np
import scipy
#####################################################################################################
# Problem 5
# Use the programs you just wrote with the same matrix and using the same settings to answer the following.
# (a) How many iterations are required for each method to reach the stopping criterion (relative error) for tol=10e-6 and tol=10e-8? Also:
#	- For each method, how does the absolute error (from the prev. question) w/ tol=10e-6 compare to the relative error?
#	- Which method required the fewest iterations?
#	- What do you observe about reaching a tighter convergence tolerance?
#####################################################################################################
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
# We redefine the solvers to run using relative tolerance instead of the absolute tolerance defined in the solvers created for problem 3.
def jacobi_solver_rel(A, b, tol=1e-6):
	n=b.size
	if(b.size!=A.shape[0] or b.size!=A.shape[1]):
		sys.exit('dimensions of A and b do not agree')
	x_old=np.transpose(np.matrix(np.zeros(n)))
	D=np.diag(np.diag(A))
	D_inv=np.linalg.inv(D)
	error=1
	counter=0
	while(error>tol):
		x_new=D_inv*(D-A)*x_old+D_inv*b
		#print(x_new)
		error=np.linalg.norm(abs(x_new-x_old)/abs(x_new))
		x_old=x_new
		counter=counter+1
	print('counter= '+str(counter))
	print('relative error = '+str(error))
	return(x_new)
#####################################################################################################
# b) Gauss-Seidel Method
def gs_solver_rel(A, b, tol=1e-6):
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
	error=1
	counter=0
	while(error>tol):
		x_new=DL_inv*(-U*x_old+b)
		#print(x_new)
		error=np.linalg.norm(abs(x_new-x_old)/abs(x_new))
		x_old=x_new
		counter=counter+1
	print('counter= '+str(counter))
	print('relative error = '+str(error))
	return(x_new)
#####################################################################################################
# c) SOR Method
def sor_solver_rel(A, b, tol=1e-6, w=1.1):
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
	DL_inv=np.linalg.inv(D+w*L)
	error=1
	counter=0
	while(error>tol):
		x_new=DL_inv*(((1-w)*D-w*U)*x_old+w*b)
		#print(x_new)
		error=np.linalg.norm(abs(x_new-x_old)/abs(x_new))
		x_old=x_new
		counter=counter+1
	print('counter= '+str(counter))
	print('relative error = '+str(error))
	return(x_new)
#####################################################################################################
# Script that executes when the program file is called from the command line.
n=input("Enter number of equations in system: ")
tol=input("Enter tolerance: ")
A=matrix_build(n)
b=b_build(n)
print('JACOBI SOLVER:')
jacobi_ans=jacobi_solver_rel(A,b,tol)
print('jacobi answer:')
print(jacobi_ans)
print('\n')

print('GAUSS-SEIDEL SOLVER:')
gs_ans=gs_solver_rel(A,b,tol)
print('gs answer:')
print(gs_ans)
print('\n')

print('SOR SOLVER:')
sor_ans=sor_solver_rel(A,b,tol)
print('sor answer:')
print(sor_ans)
print('\n')




