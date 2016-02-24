#################################################################################
# Milos Atz
# NE155 Homework 3
#################################################################################
import math
import numpy
import scipy
import matplotlib.pyplot as plt
#################################################################################
# Problem 4: Write a code that performs Composite Simpson's 3/8 rule to compute
# the integral. I is the value of the integral; a and b are the bounds; n is the
# number of points used in the integration (must be divisible by 3).

# The description of Composite Simpson is given on page 11 of the notes. It does
# Simpson's 3/8 rule on each subinterval. The endpoints of each interval are
# re-used.

def CompSimp38(a,b,n):
	if((n/3.0).is_integer()==False):
		print('please make sure n is divisible by 3, kthx')
	else:
		h=(b-a)/float(n)
		xPoints = [None]*(n+1)
		fx=[None]*(n+1)
		for i in range(0,n+1):
			xPoints[i]=a+i*h
			fx[i]=xPoints[i]/math.sqrt(xPoints[i]**2-4)
# I know this isn't the best way to do this, but I'm tired and need to get this
# done...and it works, so that's good.
		numSegments=n/3
		scndTermPts=[None]*numSegments
		for i in range(0,numSegments):
			scndTermPts[i]=fx[3*(i+1)-2]+fx[3*(i+1)-1]
		thrdTermPts=[None]*(numSegments-1)
		for i in range(0,numSegments-1):
			thrdTermPts[i]=fx[3*(i+1)]
		I=(3*h/8)*(fx[0]+3*sum(scndTermPts)+2*sum(thrdTermPts)+fx[n])
		return(I)

#################################################################################
# To determine the rate of convergence as a function of h, utilize the definition
# of "rate of convergence" from wikipedia:
# https://en.wikipedia.org/wiki/Rate_of_convergence
# We can plot this against h to determine how it behaves with decreasing h (in
# other words, increasing number of steps.

ans=2*math.sqrt(3)-math.sqrt(5)
npts=100
h=[None]*npts
compSimp=[None]*npts
for i in range(0,len(h)):
	n=3*(i+1)
	h[i]=1.0/n
	compSimp[i]=CompSimp38(3,4,n)

rateOfConv=[None]*npts
for i in range(1, len(h)):
	rateOfConv[i]=abs(compSimp[i]-ans)/abs(compSimp[i-1]-ans)

fig=plt.plot(h, rateOfConv, label='Rate of Convergence')
plt.ylabel('RoC')
plt.xlabel('h')
#plt.show()
plt.savefig('p4.png', bbox_inches='tight')

#################################################################################
# The rate of convergence is 1. The sequence is said to converge sublinearly.
#################################################################################








