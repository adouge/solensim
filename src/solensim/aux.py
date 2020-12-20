"""Units & auxilliary functions."""
import numpy as np


mm = 10**(-3)
cm = 10**(-2)
MeV = 10**6
fm = 10**(-15)


def wip():
    print("WIP feature, not implemented.")


def floor(x, unit):
    return np.floor(x/unit)*unit


def ceil(x, unit):
    return np.ceil(x/unit)*unit
