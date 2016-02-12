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
# Problem 1

pi=math.pi

def fxn(x):
	f=math.cos(pi*x/2.0)+((x**2)/2.0)
	return f

def Lagrange(data, x):
	if(len(data)==4):
		L=[None]*4
		L[0]=(x-data[1])*(x-data[2])*(x-data[3])/((data[0]-data[1])*(data[0]-data[2])*(data[0]-data[3]))
		L[1]=(x-data[0])*(x-data[2])*(x-data[3])/((data[1]-data[0])*(data[1]-data[2])*(data[1]-data[3]))
		L[2]=(x-data[0])*(x-data[1])*(x-data[3])/((data[2]-data[0])*(data[2]-data[1])*(data[2]-data[3]))
		L[3]=(x-data[0])*(x-data[1])*(x-data[2])/((data[3]-data[0])*(data[3]-data[1])*(data[3]-data[2]))
		return(L)
	else:
		print('error! not enough points in data set')

def P(x, lag_fxn, data, obj_fxn):
	if(isinstance(lag_fxn, types.FunctionType) and isinstance(obj_fxn, types.FunctionType)):
		L = lag_fxn(data, x)
		P=obj_fxn(data[0])*L[0] + obj_fxn(data[1])*L[1] + obj_fxn(data[2])*L[2] + obj_fxn(data[3])*L[3]
		return(P)

#################################################################################
# Part (c)

data=[0, 2, 3, 4]
xStart=-0.5
xEnd=4.5
nPoints=100
xIncr=(xEnd-xStart)/nPoints
xPlot=[None]*nPoints
pPlot=[None]*nPoints
fPlot=[None]*nPoints
for i in range(0,len(xPlot)):
	xPlot[i]=xStart+i*xIncr
	pPlot[i]=P(xPlot[i], Lagrange, data, fxn)
	fPlot[i]=fxn(xPlot[i])

# Gives some weird floating number errors, but it's pretty close, so we'll call it good
# print(xPlot[0:25])

figF=plt.plot(xPlot, fPlot, label='f(x)')
figP=plt.plot(xPlot, pPlot, label='P(x)')
plt.ylabel('f(x)')
plt.xlabel('x')
plt.legend(loc=4)
plt.savefig('p1c.png', bbox_inches='tight')
plt.clf()

#################################################################################
# Part (d)

data=[0, 1, 2.5, 4]
xStart=-0.5
xEnd=4.5
nPoints=100
xIncr=(xEnd-xStart)/nPoints
xPlot=[None]*nPoints
pPlot=[None]*nPoints
fPlot=[None]*nPoints
for i in range(0,len(xPlot)):
	xPlot[i]=xStart+i*xIncr
	pPlot[i]=P(xPlot[i], Lagrange, data, fxn)
	fPlot[i]=fxn(xPlot[i])

figF=plt.plot(xPlot, fPlot, label='f(x)')
figP=plt.plot(xPlot, pPlot, label='P(x)')
plt.ylabel('f(x)')
plt.xlabel('x')
plt.legend(loc=4)
plt.savefig('p1d.png', bbox_inches='tight')




