#################################################################################
# Milos Atz
# NE155 Homework 2
#################################################################################
import math
import numpy
import scipy
import types
#################################################################################
# Problem 2

pi=math.pi

def error(x, data):
	if(len(data)==4):
		error=(pi**4/384)*(x-data[0])*(x-data[1])*(x-data[2])*(x-data[3])
		return(error)

data=[0, 2, 3, 4]
xStart=0.0
xEnd=4.0
nPoints=100
xIncr=(xEnd-xStart)/nPoints
xPlot=[None]*(nPoints+1)
err=[None]*(nPoints+1)
for i in range(0,len(xPlot)):
	xPlot[i]=xStart+i*xIncr
	err[i]=abs(error(xPlot[i], data))

print(max(err))

