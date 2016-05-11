#####################################################################################################
# ANALYTICAL TEST
# This script calculates and plots the 1D flux for a plane source in a slab geometry. Under normal circumstances, this should match the output of the MC script when the inputs are set identical (assuming the tally results of the MC simulation are multiplied by the source strength).
#####################################################################################################
import numpy as np
import math
from math import exp
import random
import sys
#####################################################################################################
#import inp_test
#totalParticles = inp_test.totalParticles
#source = inp_test.sourceStrength
#xs = inp_test.xs
#rBound = inp_test.rBound
#lBound = inp_test.lBound
#xSource = inp_test.xSource
#####################################################################################################
# Function to calculate flux based on the analytical solution to the diffusion equation in 1D slab for a plane source
def flux(x, xSource, rBound, lBound, S, Lsquared, D):
	if(x > xSource):
		D = D[2]
		L = math.sqrt(Lsquared[2])
		rBound = rBound+0.71*D
		flux = (S*L/2.0/D/(1+exp(-2*(rBound-xSource)/L)))*(exp(-(x-xSource)/L)-exp(-2*(rBound-xSource)/L)*exp((x-xSource)/L))
	if(x < xSource):
		D = D[1]
		L = math.sqrt(Lsquared[1])
		lBound = lBound-0.71*D
		flux = -(S*L/2.0/D/(1+exp(-2*(lBound-xSource)/L)))*(exp(-(x-xSource)/L)-exp(-2*(lBound-xSource)/L)*exp((x-xSource)/L))
	if(x==xSource):
		print('flux undefined at x = 0')
	return(flux)
#####################################################################################################
# Return a result that can be plotted
def ana_flux(xs, xSource, rBound, lBound, sourceStrength):
	# TEST the boundaries for correctness
	if(lBound > rBound):
		sys.exit('left bound must be located to the left of the right bound')
	if(xSource < lBound or xSource > rBound):
		sys.exit('source must be located inside problem bounds')
	# Calculate absorption from the XS in the input file
	xs['a'] = {1:xs['t'][1]-xs['s'][1], 2:xs['t'][2]-xs['s'][2]}
	# Because in the MC simulation, the direction is either forward or backward, the scattering angle is always either 180 or 0 degrees. The transport xs is normally given by xs_t - 3*mu*xs_s (where mu is the average cosine of the scattering angle). In this case, mu is either -1 or 1, thus making the average mu = 0.
	# Define the diffusion coefficient in terms of XS
	D = {}
	D[1] = 1.0/(3*(xs['t'][1]))
	D[2] = 1.0/(3*(xs['t'][2]))
	# Define the diffusion length
	Lsq = {}
	Lsq[1] = D[1]/xs['a'][1]
	Lsq[2] = D[2]/xs['a'][2]
	# With the constants defined, we can evaluate the analytical solution of the flux over the desired space.
	nPoints = int((rBound-lBound)*100)	# increments of 0.01 for solution
	xPoints = [0]*(nPoints+1)
	xPoints[0] = lBound
	for i in range(1, nPoints+1):
		xPoints[i]=round(xPoints[i-1]+0.01, 3)
	xPoints.remove(xSource)
	fluxSoln = [0]*(nPoints)
	fluxSoln[0] = flux(xPoints[0], xSource, rBound, lBound, sourceStrength, Lsq, D)
	for i in range(0, nPoints):
		fluxSoln[i] = flux(xPoints[i], xSource, rBound, lBound, sourceStrength, Lsq, D)
	analyticalResult = [xPoints, fluxSoln]
	return(analyticalResult)
##########################################################################
# plot the solutions for (c) and (d)
#import matplotlib.pyplot as plt
#fig=plt.plot(xPoints, fluxSoln)
#plt.ylabel('flux')
#plt.xlabel('x')
#plt.show()
##plt.savefig('p1.png', bbox_inches='tight')
#plt.clf()
