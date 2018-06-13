from __future__ import division


def anint(double val):
    if val >= 0.5:
        return 1.0
    elif val <= -0.5:
        return -1.0
    else:
        return 0.0
