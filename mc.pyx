from __future__ import division

from copy import deepcopy
from math import exp
import random


cdef int kb = 1

class MonteCarlo():
    """
    Parameters
    ----------
    system : System
        System on which to perform MC
    temperature : float
        Temperature of the MC simulation
    dx : float
        Maximum displacement. Will be adjusted to achieve a target acceptance
        probability.
    target : float
        Target acceptance probability
    seed : int, optional, default=12345
        Seed for random number generator
    """
    def __init__(self, system, double temperature, double dx, double target,
                 int seed=12345):
        self.system = system
        self.dx = dx
        self.target = target

        self.temperature = temperature
        self.beta = 1 / (kb * temperature)

        random.seed(seed)

        self.n_accept = 0

    def take_step(self):
        cdef int id
        cdef double pe_old, pe_new, delta_pe, rand_value

        for id in range(self.system.n):
            pe_old = self.system.calc_energy(id=id)
            xyz_old = deepcopy(self.system.traj.xyz[0])
            self.system.displace_coords(id, self.dx)
            pe_new = self.system.calc_energy(id=id)

            if pe_new:
                delta_pe = pe_new - pe_old
                if delta_pe < 0:
                    self.n_accept += 1
                else:
                    rand_value = random.random()
                    if exp(-self.beta * delta_pe) > rand_value:
                        self.n_accept += 1
                    else:
                        self.system.traj.xyz[0] = xyz_old
            else:
                self.system.traj.xyz[0] = xyz_old

    def relax(self, int steps, int adjust_freq=100):
        """
        Parameters
        ----------
        adjust_freq : int
            How often to adjust `dx` to achieve the target acceptance probability in
            number of steps.
        """
        cdef int i
        cdef double prob, potential

        for i in range(steps):
            self.take_step()
            self.system.check_nlist()

            if i != 0 and i % adjust_freq == 0:
                prob = self.n_accept / (adjust_freq * self.system.n)
                if prob < self.target:
                    self.dx /= 1.025
                elif prob > self.target:
                    self.dx *= 1.025

                if self.dx > self.system.skin / 2:
                    self.dx = self.system.skin / 2

                potential = self.system.calc_energy(id=0, total=True)

                print('Relax: {} of {}\tdx: {:.6f}\tprob: {:.6f}\ttarget_prob: {}'
                      '  PE: {:.5f}'.format(i, steps, self.dx, prob, self.target,
                                            potential / self.system.n))

                self.n_accept = 0

    def run(self, str traj_filename, int steps, int freq=1000):
        cdef int i
        cdef double prob

        for i in range(steps):
            self.take_step()
            self.system.check_nlist()

            if i != 0 and i % freq == 0:
                prob = self.n_accept / (freq * self.system.n)

                potential = self.system.calc_energy(id=0, total=True)

                print('Run: {} of {}\tPE: {:.5f}\tprob: {:.6f}'
                      ''.format(i, steps, potential / self.system.n, prob))

                self.n_accept = 0
                self.system.write_xyz(traj_filename)
