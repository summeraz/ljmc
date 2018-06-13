from __future__ import division

from utils import *


### Define Lennard-Jones Parameters ###
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

### Define System Parameters ###
n_particles = 256
number_density = 0.7

### Define Monte Carlo Parameters ###
temperature = 1.2   # Temperature of the simulation
dx = 0.1            # Initial maximum displacement
target = 0.5        # Target acceptance probabality
n_relax = 10000     # Number of timesteps to relax from initial configuration
n_mc = 0       # Total number of MC steps

#############################################################################

#######
# RUN #
#######

# Create the force field
forcefield = ForceField(sigma=sigma, epsilon=epsilon, cutoff=cutoff)

# Create the system
system = System(n_particles, number_density, forcefield)

# Initialize the neighborlist
system.build_nlist(skin=0.5)

# Create Monte Carlo instance
mc = MonteCarlo(system=system, dx=dx, temperature=temperature, target=target)

'''
import cProfile, pstats
pr = cProfile.Profile()
pr.enable()
'''

# Relax the system (and optimize `dx`)
mc.relax(n_relax, adjust_freq=20)

'''
pr.disable()
f = open('run-mdtraj.prof', 'a')
sortby = 'cumulative'
pstats.Stats(pr, stream=f).strip_dirs().sort_stats(sortby).print_stats()
f.close()
'''

# Monte Carlo production run
#mc.run(n_mc)
