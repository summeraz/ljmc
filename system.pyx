from __future__ import division

from copy import deepcopy
import itertools
import os
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

        box_length = (n / density) ** (1 / 3)
        box = mb.Box(lengths=np.ones(3) * box_length)
        self.box = box
        self.box_length = box_length

        particles_per_len = round(n ** (1 / 3))
        pattern = mb.Grid3DPattern(particles_per_len, particles_per_len,
                                   particles_per_len)
        pattern.scale(box_length)
        compound = mb.Compound()
        for pos in pattern:
            compound.add(mb.Particle(name='_LJ', pos=pos))
        compound.periodicity = box.lengths
        self.traj = compound.to_trajectory()

        self.print_params()

    def calc_energy(self, int id, double inner_cutoff=0.65, total=False):
        """
        Parameters
        ----------
        id : int
            Return PE of the particle with this ID (if total is `False`)
        inner_cutoff : float, optional, default=0.65
            Close-range cutoff for inter-particle distances to avoid overflow
        total : bool, optional, default=False
            Return the total PE of the system
        """
        cdef int i
        cdef double cutoff, epsilon, sigma, dist, pe

        cutoff = self.forcefield.cutoff
        epsilon = self.forcefield.epsilon
        sigma = self.forcefield.sigma

        pe = 0
        for i in range(self.n):
            if total or id == i:
                atom_pairs = list(itertools.product([i], self.nlist[i]))
                dists = md.compute_distances(self.traj, atom_pairs)[0]
                if not total and min(dists) < inner_cutoff:
                    return False
                else:
                    for dist in dists:
                        pe += 4 * epsilon * ((sigma / dist)**12 - (sigma / dist)**6)
        if total:
            pe /= 2

        return pe

    def displace_coords(self, int id, double dx):
        """
        Parameters
        ----------
        id : int
            ID of particle to displace
        dx : float
            Maximum displacement distance in each dimension
        """
        cdef int k
        cdef double dx_temp, new_pos

        for k in range(3):
            dx_temp = dx * (2.0 * random.random() - 1.0)
            new_pos = self.traj.xyz[0, id, k] + dx_temp
            new_pos -= self.box_length * anint(new_pos * (1 / self.box_length))
            self.traj.xyz[0, id, k] = new_pos

    def build_nlist(self, double skin):
        """
        Parameters
        ----------
        skin : float
            Skin distance for the neighborlist
        """
        cdef double total_range

        self.skin = skin
        self.xyz_old = deepcopy(self.traj.xyz[0])

        total_range = skin + self.forcefield.cutoff
        nlist = [md.compute_neighbors(self.traj, total_range, [id])[0]
                 for id in range(self.n)]
        self.nlist = nlist

    def check_nlist(self):
        """
        """
        cdef int k
        cdef double half_skin

        half_skin = self.skin / 2
        rebuild = False

        for xyz_new, xyz_old in zip(self.traj.xyz[0], self.xyz_old):
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
        self.traj.save_xyz('tmp.xyz')
        os.system('cat tmp.xyz >> {}'.format(filename))
        os.remove('tmp.xyz')
