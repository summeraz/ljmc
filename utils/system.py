from __future__ import division

from copy import deepcopy
import itertools
import random

import mbuild as mb
import numpy as np

from .utils import *


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
    def __init__(self, n, density, forcefield):
        self.n = n
        self.density = density
        self.forcefield = forcefield

        box_length = (n / density) ** (1 / 3)
        box = mb.Box(lengths=np.ones(3) * box_length)
        self.box = box
        self.box_length = box_length

        particle = mb.Particle(name='_LJ')
        packed_box = mb.fill_box(particle, n_compounds=n, box=box, overlap=0.9,
                                 edge=0.9)
        self.xyz = packed_box.xyz

        self.print_params()

    def calc_energy(self, inner_cutoff=0.75, id=None):
        """
        Parameters
        ----------
        inner_cutoff : float, optional, default=0.75
            Close-range cutoff for inter-particle distances to avoid overflow
        id : int, optional, default=None
            Only return PE of the particle with this ID
        """
        cutoff = self.forcefield.cutoff
        epsilon = self.forcefield.epsilon
        sigma = self.forcefield.sigma

        pe = 0
        for i, xyz1 in enumerate(self.xyz):
            if not id or id == i:
                for j in self.nlist[i]:
                    xyz2 = self.xyz[j]
                    rij = []
                    for k in range(3):
                        rij.append(xyz1[k] - xyz2[k])
                        pbc = self.box_length * anint(rij[k] * (1 / self.box_length))
                        rij[k] -= pbc
                    r = np.linalg.norm(rij)
                    if r < inner_cutoff:
                        return False
                    elif r < cutoff:
                        pe += 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
        if not id:
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
        for k in range(2):
            dx_temp = dx * (2.0 * random.random() - 1.0)
            new_pos = self.xyz[id, k] + dx_temp
            new_pos -= self.box_length * anint(new_pos * (1 / self.box_length))
            self.xyz[id, k] = new_pos

    def build_nlist(self, skin):
        """
        Parameters
        ----------
        skin : float
            Skin distance for the neighborlist
        """
        self.skin = skin
        self.xyz_old = deepcopy(self.xyz)
        nlist = [[] for _ in self.xyz]
        total_range = skin + self.forcefield.cutoff

        for i, xyz in enumerate(self.xyz):
            for j, xyz2 in enumerate(self.xyz):
                if i != j:
                    rij = []
                    for k in range(3):
                        rij.append(xyz[k] - xyz2[k])
                        pbc = self.box_length * anint(rij[k] * (1 / self.box_length))
                        rij[k] -= pbc
                    if np.linalg.norm(rij) < total_range:
                        nlist[i].append(j)

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
