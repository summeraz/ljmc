from __future__ import division


class ForceField():
    def __init__(self, sigma, epsilon, cutoff):
        self.sigma = sigma
        self.epsilon = epsilon
        self.cutoff = cutoff
