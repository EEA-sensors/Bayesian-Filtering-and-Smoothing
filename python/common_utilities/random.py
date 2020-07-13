"""
A random utility file used for reproducibility across languages (Octave, Julia, etc). Enforces the use of
Mersen-Twister 1997 and of some algorithms.
"""
import numpy as np
import numpy.random as rd
from scipy.special import erfcinv

__all__ = ["RandomState"]

SQRT_2 = 2 ** 0.5


class RandomState(rd.RandomState):
    """Enforces MT1997 and reproducibility with Matlab"""

    def __init__(self, seed=0):
        if not isinstance(seed, int):
            raise TypeError(f"seed must be an int, '{seed}' of type '{type(seed)}' was passed")
        super(RandomState, self).__init__(seed)

    def randn(self, *args):
        if args:
            uniforms = self.rand(*args[::-1]).T
        else:
            uniforms = self.rand()
        return SQRT_2 * erfcinv(2 * uniforms)

    def choice(self, n, k, p=None):
        if p is None:
            p = np.full(n, 1/n)
        cs = np.cumsum(p)
        cs = cs / cs[-1]
        uniforms = self.rand(k)
        return np.searchsorted(cs, uniforms, "left")

