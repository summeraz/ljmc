from __future__ import division


def anint(double val):
    if val >= 0.5:
        return 1
    elif val <= -0.5:
        return -1
    else:
        return 0
