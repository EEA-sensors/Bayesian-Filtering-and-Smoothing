"""
A random utility file used for reproducibility across languages (Octave, Julia, etc). Enforces the use of
Mersen-Twister 1997
"""
import numpy.random as rd
from scipy.special import erfinv

__all__ = ["RandomState"]

SQRT_2 = 2 ** 0.5


class RandomState(rd.RandomState):
    """Enforces MT1997"""

    def __init__(self, seed=0):
        if not isinstance(seed, int):
            raise TypeError(f"seed must be an int, '{seed}' of type '{type(seed)}' was passed")
        super(RandomState, self).__init__(seed)

    def randn(self, *args):
        uniforms = self.rand(*args)
        return SQRT_2 * erfinv(2 * uniforms - 1)
