from __future__ import division

import pyximport; pyximport.install()

from forcefield import *
from mc import *
from system import *
from utils import *


### Define Lennard-Jones Parameters ###
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

### Define System Parameters ###
n_particles = 125
number_density = 0.5

### Define Monte Carlo Parameters ###
temperature = 1.2   # Temperature of the simulation
dx = 0.1            # Initial maximum displacement
target = 0.5        # Target acceptance probabality
n_relax = 2500      # Number of timesteps to relax from initial configuration
n_mc = 25000        # Total number of MC steps

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

# Relax the system and optimize `dx`
mc.relax(n_relax, adjust_freq=50)

# Monte Carlo production run
mc.run(traj_filename='traj.xyz', steps=n_mc, freq=100)
