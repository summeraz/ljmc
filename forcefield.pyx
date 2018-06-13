from __future__ import division


class ForceField():
    def __init__(self, double sigma, double epsilon, double cutoff):
        self.sigma = sigma
        self.epsilon = epsilon
        self.cutoff = cutoff
