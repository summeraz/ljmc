from __future__ import division

from copy import deepcopy
import itertools
import random

import mbuild as mb
import mdtraj as md
import numpy as np

from utils import *


class System():
    """
    Parameters
    ----------
    n : int
        Number of particles in the system
    density : float
        Number density of particles in the system
    forcefield : ForceField
        Force field defining particle interactions
    """
    def __init__(self, int n, double density, forcefield):
        self.n = n
        self.density = density
        self.forcefield = forcefield

        cdef float box_length = (n / density) ** (1 / 3)
        box = mb.Box(lengths=np.ones(3) * box_length)
        self.box = box
        self.box_length = box_length

        particle = mb.Particle(name='_LJ')
        packed_box = mb.fill_box(particle, n_compounds=n, box=box, overlap=0.9,
                                 edge=0.9)
        self.traj = packed_box.to_trajectory()
        self.xyz = self.traj.xyz[0]

        self.print_params()

    def calc_energy(self, int id, double inner_cutoff=0.75, total=False):
        """
        Parameters
        ----------
        id : int
            Return PE of the particle with this ID (if total is `False`)
        inner_cutoff : float, optional, default=0.75
            Close-range cutoff for inter-particle distances to avoid overflow
        total : bool, optional, default=False
            Return the total PE of the system
        """
        cdef double cutoff, epsilon, sigma
        cutoff = self.forcefield.cutoff
        epsilon = self.forcefield.epsilon
        sigma = self.forcefield.sigma

        cdef int i
        cdef double dist, pe
        pe = 0
        for i in range(self.n):
            if total or id == i:
                atom_pairs = list(itertools.product([i], self.nlist[i]))
                dists = md.compute_distances(self.traj, atom_pairs)[0]
                if min(dists) < inner_cutoff:
                    return False
                else:
                    for dist in dists:
                        pe += 4 * epsilon * ((sigma / dist)**12 - (sigma / dist)**6)
        if total:
            pe /= 2

        return pe

    def displace_coords(self, id, dx):
        """
        Parameters
        ----------
        id : int
            ID of particle to displace
        dx : float
            Maximum displacement distance in each dimension
        """
        for k in range(3):
            dx_temp = dx * (2.0 * random.random() - 1.0)
            new_pos = self.xyz[id, k] + dx_temp
            new_pos -= self.box_length * anint(new_pos * (1 / self.box_length))
            self.traj.xyz[0, id, k] = new_pos

    def build_nlist(self, skin):
        """
        Parameters
        ----------
        skin : float
            Skin distance for the neighborlist
        """
        self.skin = skin
        self.xyz_old = deepcopy(self.xyz)

        total_range = skin + self.forcefield.cutoff
        nlist = [md.compute_neighbors(self.traj, total_range, [id])[0]
                 for id in range(self.n)]
        self.nlist = nlist

    def check_nlist(self):
        """
        """
        half_skin = self.skin / 2
        rebuild = False

        for xyz_new, xyz_old in zip(self.xyz, self.xyz_old):
            rij = []
            for k in range(3):
                rij.append(xyz_new[k] - xyz_old[k])
                pbc = self.box_length * anint(rij[k] * (1 / self.box_length))
                rij[k] -= pbc
            if np.linalg.norm(rij) > half_skin:
                rebuild = True
                break

        if rebuild:
            self.build_nlist(skin=self.skin)

    def print_params(self):
        """
        """
        print('=======================')
        print('   System Parameters')
        print('=======================')
        print('n_particles = {}'.format(self.n))
        print('number density = {}'.format(self.density))
        print('box dimensions = {:.2f} x {:.2f} x {:.2f}'.format(*self.box.lengths))
        print('=======================')

    def write_xyz(self, filename):
        """
        Parameters
        ----------
        filename : str
            Filename to write coordinates to
        """
        pass
