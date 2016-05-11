#####################################################################################################
# TEST RUN FILE
# This contains the input data and specifications for the test simulation against the analytical solution of the diffusion equation for plane source with vacuum BC
#####################################################################################################
# INPUT DATA
#####################################################################################################
# Enter the number of particles to simulate (recommend 1e4 - 1e5)
totalParticles = 1e4
# The source strength; used in the analytical solution and to scale the MC solution
sourceStrength = 10
# Input cross sections; s: scattering, t: total
xs = {}
xs['s'] = {1:0.5, 2:0.5}
xs['t'] = {1:1.0, 2:1.0}
# Specify the location of the left and right bounds, as well as the location of the plane source. Make sure the plane source is located within the bounds and that the left bound is less than the right bound.
rBound = 6
lBound = -4
xSource = 1
xInterface = 0
# Specify the width of the bins for the tally mesh (recommend 0.1)
stepSize = 0.25
# Include analytical solution to the diffusion equation?
incl = 'yes'
#####################################################################################################
# END OF INPUT
#####################################################################################################
import matplotlib.pyplot as plt
if(incl == 'yes'):
	# Force the source point and interface to be the same location
	xInterface = xSource
	# Import the monte carlo solver; use to solve for the flux in the space defined by the input parameters. Note that the interface between the two regions of the problem MUST be defined as the position of the source for the analytical solution to make sense in conjunction with the mc solution
	import mc
	mcSoln = mc.mc(totalParticles, xs, rBound, lBound, xInterface, xSource, sourceStrength, stepSize)
	print('maximum relative error = '+str(max(mcSoln[3])))
	# Import the analytical solution script and functions; use to solve for the flux in the space defined by the input parameters
	import analytical_test
	analyticalSoln = analytical_test.ana_flux(xs, xSource, rBound, lBound, sourceStrength)
	fig_an=plt.plot(analyticalSoln[0], analyticalSoln[1], label='analytical')
#####################################################################################################
else:
	# Import the monte carlo solver; use to solve for the flux in the space defined by the input parameters. Note that the interface between the two regions of the problem MUST be defined as the position of the source for the analytical solution to make sense in conjunction with the mc solution
	import mc
	mcSoln = mc.mc(totalParticles, xs, rBound, lBound, xInterface, xSource, sourceStrength, stepSize)
	print('maximum relative error = '+str(max(mcSoln[3])))
# Compare the analytical and monte carlo results by plotting them together in space.
fig_mc=plt.plot(mcSoln[0], mcSoln[1], 'ro', label='mc')
fig_mc=plt.errorbar(mcSoln[0], mcSoln[1], yerr=mcSoln[3], linestyle="None")
#####################################################################################################
plt.ylabel('flux')
plt.xlabel('x')
plt.legend()
plt.figtext(0.14, 0.86, 'Source Location = '+str(xSource))
plt.figtext(0.14, 0.82, 'Interface Location = '+str(xInterface))
plt.figtext(0.14, 0.78, 'Number of histories = '+str(int(totalParticles)))
plt.figtext(0.14, 0.74, 'Maximum relative error = '+str(round(max(mcSoln[3]),4)))
plt.show()
#plt.savefig('result.png', bbox_inches='tight')
plt.clf()
#####################################################################################################
# The diffusion and MC result show good agreement away from the source. The disagreement at the source is due to two things:
# First, in the analytical solution, the point corresponding to the source point is removed from the set of points on which the flux is calculated (because it is not defined at the source). This creates the jagged line sometimes seen in the analytical solution.
# In addition, the reason the analytical solution might not match up perfectly with the monte carlo result is due to the limitations of diffusion theory, which assumes that the flux is slowly varying in space and mostly contributed to by scattering. Near sources, both of these assumptions break down. Further, the analytical solution neglects particles crossing the boundary (which is also the position of the source) by assuming that the current approaching the boundary trends to S/2.
# Also, the mc solution requires binning the collisions within the bin interval; if not many collisions happen, there can be some disagreement between the recorded tally and the analytical result. To improve precision, can use coarser bins.
#####################################################################################################
