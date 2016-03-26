#####################################################################################################
# Milos Atz
# NE155 Homework 5
#####################################################################################################
import math
import numpy as np
import scipy
#####################################################################################################
# Problem 4
# Write a program to implement the following iterative methods for a matrix with n unknowns.
# (a) Jacobi method
# (b) Gauss Seidel method
# (c) SOR method
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
# a) Jacobi method
# Strategy: use a while loop to iterate while some convergence equation between x_old and x_new is greater than the input tolerance. The initial values for x are defined as x_old. In every while loop, x_new is calculated based on x_old. The error is then calculated. If the error tolerance is met, the while loop ends; if not, x_old = x_new and the loop repeats. I should implement an iteration counter to count the number of loops.
def jacobi_solver(A, b, tol=1e-6):
	n=b.size
	if(b.size!=A.shape[0] or b.size!=A.shape[1]):
		sys.exit('dimensions of A and b do not agree')
	x_old=np.transpose(np.matrix(np.zeros(n)))
	D=np.diag(np.diag(A))
	D_inv=np.linalg.inv(D)
	conv=1
	counter=0
	while(conv>tol):
		x_new=D_inv*(D-A)*x_old+D_inv*b
		#print(x_new)
		conv=np.linalg.norm(x_new-x_old)
		x_old=x_new
		counter=counter+1
	print('counter= '+str(counter))
	print('absolute error = '+str(conv))
	return(x_new)
#####################################################################################################
# b) Gauss-Seidel Method
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
		#print(x_new)
		conv=np.linalg.norm(x_new-x_old)
		x_old=x_new
		counter=counter+1
	print('counter= '+str(counter))
	print('absolute error = '+str(conv))
	return(x_new)
#####################################################################################################
# c) SOR Method
def sor_solver(A, b, w=1.1, tol=1e-6):
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
	conv=1
	counter=0
	while(conv>tol):
		x_new=DL_inv*(((1-w)*D-w*U)*x_old+w*b)
		#print(x_new)
		conv=np.linalg.norm(x_new-x_old)
		x_old=x_new
		counter=counter+1
	print('counter= '+str(counter))
	print('absolute error = '+str(conv))
	return(x_new)
#####################################################################################################
# Script that executes when the program file is called from the command line.
n=input("Enter number of equations in system: ")
A=matrix_build(n)
b=b_build(n)
print('JACOBI SOLVER:')
jacobi_ans=jacobi_solver(A,b)
print('jacobi answer:')
print(jacobi_ans)
print('\n')

print('GAUSS-SEIDEL SOLVER:')
gs_ans=gs_solver(A,b)
print('gs answer:')
print(gs_ans)
print('\n')


print('SOR SOLVER:')
sor_ans=sor_solver(A,b)
print('sor answer:')
print(sor_ans)
print('\n')
