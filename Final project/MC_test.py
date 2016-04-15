#####################################################################################################
import numpy as np
import math
import random
#####################################################################################################
# Sample particle direction - this returns just an integer letting you know whether the particle is moving in the positive or negative direction.
def sample_direction():
	if(random.random() >= 0.5):
		direction = 1
	else:
		direction = -1
	return(direction)
#####################################################################################################
# Sample path length for the region that the neutron is in. If the neutron is on the boundary, sample the path length for the region that the neutron is entering.
def sample_path_length(x, dir, xs):
	if(x == 0):
		x = x+dir		# Probably can do this one in terms of region as well!
	if(x < 0):
		xs_t = xs['t'][1]
		xs_s = xs['s'][1]
	if(x > 0):
		xs_t = xs['t'][2]
		xs_s = xs['s'][2]
	pathLength = -math.log(random.random())/xs_t
	return(pathLength)
#####################################################################################################
# Sample collision
def sample_collision(x, xs):
	if(x < 0):
		xs_t = xs['t'][1]
		xs_s = xs['s'][1]
	if(x > 0):
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
## Calculate what needs to be added to the running tally according to the rules of implicit capture.
#def add_tally(tally, region, xs, weight):
#	xs_t = xs['t'][region]
#	xs_s = xs['s'][region]
#	newTally = weight*(xs_t-xs_s)/xs_t
#	totTally = tally[region]+newTally
#	return(totTally)
#####################################################################################################
## Reduce the weight of the particle with every collision according to the rules of implicit capture.
#def reduce_weight(weight, xs, region):
#	xs_t = xs['t'][region]
#	xs_s = xs['s'][region]
#	newWeight = weight*(1-((xs_t-xs_s)/xs_t))
#	return(newWeight)
#####################################################################################################
# Define XS and weight data
xs = {}
xs['s'] = {1:0.5, 2:0.75}
xs['t'] = {1:1.0, 2:0.90}
#w = {}
#w['nom'] = {1:1.0, 2:2.0}
#w['max'] = {1:2.5, 2:5.0}
#w['min'] = {1:0.4, 2:0.8}
#####################################################################################################
totalParticles = input("Enter number of particles to simulate: ")
stepSize = input("Enter spatial bin size: ")
n = 1
#####################################################################################################
rBound = 4.0
lBound = -4.0
#####################################################################################################
# The collision tally is performed over two the two regions. They are stored separately in a dictionary. The tally that gets updated depends on the region where the collision takes place.
#tally = {1:0.0, 2:0.0}
colTally = []
colTallyA = []
colTallyS = []
#####################################################################################################
# Write the script that will run the code; this will include the summations of the tallies. The code will prompt the user for the number of particles to run.
log=open('./mc_logfile.out','w')
z = 'MC Logfile - Output'+'\n'
log.writelines(z)
while(n <= totalParticles):
	# Particle is born at x = -1 cm
	#log.writelines('particle '+str(n)+'\n')
	location = -1.0
	weight = 1.0
	while(weight == 1.0):
		dir = sample_direction()
		pathLength = sample_path_length(location,dir,xs)
		if(np.sign(location)!=np.sign(location+dir*pathLength)):
			# Cross boundary, will switch regions -> need to resample direction
			# what happens if path length = distance to boundary?
			location = 0.0
			pathLength = sample_path_length(location,dir,xs)
			#log.writelines('resampling path length to cross boundary'+'\n')
		location = location+dir*pathLength
		if(location > rBound or location < lBound):
			# particle leaked, set weight to zero
			weight = 0.0
		else:
			collision=sample_collision(location, xs)
			#log.writelines('location = '+str(location)+'\n')
			#log.writelines('collision type = '+ str(collision)+'\n')
			if(collision == 'a'):
				# tally collision, particle is terminated
				colTally.append(location)
				colTallyA.append(location)
				weight = 0.0
			else:
				# tally collision, particle is scattered: resample direction and pathLength
				colTally.append(location)
				colTallyS.append(location)
	n=n+1
#print(math.ceil(max(colTally)))
#print(math.floor(min(colTally)))
import matplotlib.pyplot as plt
n, bins, patches = plt.hist(colTally, np.linspace(math.floor(min(colTally)), math.ceil(max(colTally)), num=math.ceil(max(colTally))-math.floor(min(colTally))+1), label = 'collisions')
plt.legend()
plt.show()
plt.clf
#####################################################################################################
#n, bins, patches = plt.hist(colTallyA, np.linspace(math.floor(min(colTallyA)), math.ceil(max(colTallyA)), num=math.ceil(max(colTallyA))-math.floor(min(colTallyA))+1), label = 'absorptions')
#plt.legend()
#plt.show()
#plt.clf
#n, bins, patches = plt.hist(colTallyS, np.linspace(math.floor(min(colTallyS)), math.ceil(max(colTallyS)), num=math.ceil(max(colTallyS))-math.floor(min(colTallyS))+1), label = 'scatters')
#plt.legend()
#plt.show()
#plt.clf
#####################################################################################################
bins = [0]*(int((rBound-lBound)/stepSize)+1)
for i in range(0, len(bins)):
	bins[i]=lBound+stepSize*i
test=np.histogram(colTally, bins)
flux = [x/float(totalParticles) for x in test[0]]
flux.insert(0, 0.0)
flux.append(0.0)
points = [0]*int((rBound-lBound)/stepSize)
for i in range(0, len(points)):
	points[i]=bins[i]+0.5*stepSize
points.insert(0, lBound)
points.append(rBound)
plt.plot(points, flux)
plt.show()


