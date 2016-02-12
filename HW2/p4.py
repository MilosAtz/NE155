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
# Problem 4
# d) Make a log-log plot that displays both the input data and the function E=kh^p.


h=[5e-02, 2.5e-02, 1.25e-02, 6.25e-03, 3.125e-03, 1.5625e-03, 7.8125e-04, 3.90625e-04]
e=[1.036126e-01, 3.333834e-02, 1.375409e-02, 4.177237e-03, 1.103962e-03, 2.824698e-04, 7.185644e-05, 1.813937e-05]

p=1.790291394
k=28.43277105

xStart=min(h)
xEnd=max(h)
xIncr=min(h)
nPoints=int((xEnd-xStart)/xIncr)
xPlot=[None]*(nPoints+1)
ePlot=[None]*(nPoints+1)
for i in range(0,len(xPlot)):
	xPlot[i]=xStart+i*xIncr
	ePlot[i]=k*xPlot[i]**p

figData=plt.loglog(h, e, 'ro', label='Data points')
figEFxn=plt.loglog(xPlot, ePlot, label='Least squares approximation')
plt.ylabel('Error, f(h)')
plt.xlabel('h')
plt.legend(loc=4)
plt.savefig('p4d.png', bbox_inches='tight')
plt.clf()
