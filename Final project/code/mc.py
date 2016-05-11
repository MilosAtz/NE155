#####################################################################################################
import numpy as np
import math
import random
import sys
from bisect import bisect_left
import matplotlib.pyplot as plt
from operator import add
#####################################################################################################
# Sample particle direction - this returns just an integer letting you know whether the particle is moving in the positive or negative direction.
def sample_direction():
	if(random.random() >= 0.5):
		direction = 1
	else:
		direction = -1
	return(direction)
#####################################################################################################
# Sample the cosine of the scattering angle. At this point, it is assumed that the scattering angle is uniformly distributed for all media!
# Typically, this would range from -1 to 1. However, because direction is determined using another function, all this value must do is scale the path length by the cosine of the absolute value of the scattering angle (between 0 and 1).
def sample_angle():
	cos = random.random()
	return(cos)
#####################################################################################################
# Sample path length for the region that the neutron is in. If the neutron is on the boundary, sample the path length for the region that the neutron is entering.
def sample_path_length(x, dir, xs, interface):
	if(x == interface):
		x = x+dir		# Probably can do this one in terms of region as well!
	if(x < interface):
		xs_t = xs['t'][1]
		xs_s = xs['s'][1]
	if(x > interface):
		xs_t = xs['t'][2]
		xs_s = xs['s'][2]
	pathLength = -math.log(random.random())/xs_t
	return(pathLength)
#####################################################################################################
# Sample collision
def sample_collision(x, xs, interface):
	if(x < interface):
		xs_t = xs['t'][1]
		xs_s = xs['s'][1]
	if(x > interface):
		xs_t = xs['t'][2]
		xs_s = xs['s'][2]
	P_s = xs_s/xs_t
	P_a = (xs_t-xs_s)/xs_t
	if(random.random() >= P_s):
		collision = 'a'
	else:
		collision = 's'
	return(collision)
#####################################################################################################
# This function finds the closest number in a list to a float and returns it. If the float is equidistant from to numbers, it'll return the smaller of the two. This is used in tallying; the location of every collision is compared to the list of bins and based on that, a specific entry in the tally mesh is incremented by the weight.
def takeClosest(myList, myNumber):
	pos = bisect_left(myList, myNumber)
	if(pos == 0):
		return(myList[0])
	if(pos == len(myList)):
		return(myList[-1])
	before = myList[pos - 1]
	after = myList[pos]
	if(after - myNumber < myNumber - before):
		return(after)
	else:
		return(before)
#####################################################################################################
# This function adds the contribution of a collision in a specific bin of the tally for a single history.
def tally(ctally, x, w, bins):
	# Identify the bin in which to increment
	bin = takeClosest(bins, x)
	idx = bins.index(bin)
	# Increment the score
	ctally[idx] = ctally[idx]+w
	return(ctally)
#####################################################################################################
# Function that will run the code and collect the tallies (this is what is called in the run file)
def mc(totalParticles, xs, rBound, lBound, interface, xSource, srceStrngth, binWidth, nonanalog='off', w = 'NA'):
	# Test the boundaries for correctness
	if(lBound > rBound):
		sys.exit('left bound must be located to the left of the right bound')
	if(xSource < lBound or xSource > rBound):
		sys.exit('source must be located inside spatial bounds')
	if(interface < lBound or interface > rBound):
		sys.exit('interface between regions must be located inside spatial bounds')
	# Set up the bins for the tally mesh
	bins = [0]*(int((rBound-lBound)/binWidth)+1)
	for i in range(0, len(bins)):
		bins[i]=lBound+binWidth*i
#####################################################################################################
	# After each history, the score for each bin is squared and added to the appropriate entry in the sqTally list; this is used to calculate standard deviation.
	sqTally = [0]*len(bins)
	# After each history, the score for each bin is added to the appropriate tally in the cumTally list; this is used to calculate standard deviation and is used to calculate flux.
	cumTally = [0]*len(bins)
#####################################################################################################
	# The function offers the option to use variance reduction. If analog simulation is desired, the following will run, returning the collision tally.
	n = 1
	if(nonanalog == 'off'):
		while(n <= totalParticles):
			# The history collision tally keeps track of the location of collisions for each history. It is reset at the beginning of every history, so it is defined in the loop.
			histTally = [0]*len(bins)
			# Particle is born at x = -1 cm
			location = xSource
			weight = 1.0
			while(weight == 1.0):
				dir = sample_direction()
				pathLength = sample_path_length(location, dir, xs, interface)
				angle = sample_angle()
				if(np.sign(interface-location)!=np.sign(interface-(location+dir*pathLength*angle))):
					# Cross boundary, will switch regions -> need to resample direction
					location = interface
					pathLength = sample_path_length(location, dir, xs, interface)
				location = location+dir*pathLength*angle
				if(location > rBound or location < lBound):
					# particle leaked, set weight to zero
					weight = 0.0
				else:
					collision=sample_collision(location, xs, xSource)
					if(collision == 'a'):
						# tally collision, particle is terminated
						histTally = tally(histTally, location, weight, bins)
						weight = 0.0
					else:
						# tally collision, particle is scattered: resample direction and pathLength
						histTally = tally(histTally, location, weight, bins)
			# Increment the cumulative tally and cumulative square tallies based on the history of particle i
			for i in range(0, len(bins)):
				sqTally[i] = sqTally[i]+histTally[i]**2
				cumTally[i] = cumTally[i]+histTally[i]
			n=n+1
# Calculate the sample mean and variance using the collected tallies
	sampleMean = [0]*len(bins)
	sampleVar = [0]*len(bins)
	stdDev = [0]*len(bins)
	relErr = [0]*len(bins)
	mcFlux = [0]*len(bins)
	for i in range(0, len(bins)):
		sampleMean[i] = cumTally[i]/float(totalParticles)/binWidth
		sampleVar[i] = (sqTally[i]-2*sampleMean[i]*cumTally[i]+(sampleMean[i]**2)*totalParticles)/(totalParticles-1)
		stdDev[i] = srceStrngth*math.sqrt(sampleVar[i]/totalParticles)
		mcFlux[i] = srceStrngth*sampleMean[i]
		if(mcFlux[i]!=0):
			relErr[i] = stdDev[i]/mcFlux[i]
	return(bins, mcFlux, stdDev, relErr)
#####################################################################################################
# This final list is called back into the run file to be plotted
#####################################################################################################



