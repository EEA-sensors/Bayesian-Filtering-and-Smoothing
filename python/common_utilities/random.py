"""
A random utility file used for reproducibility across languages (Octave, Julia, etc). Enforces the use of
Mersen-Twister 1997
"""
import numpy.random as rd

__all__ = ["RandomState"]


class RandomState(rd.RandomState):
    """Enforces MT1997"""

    def __init__(self, seed=0):
        if not isinstance(seed, int):
            raise TypeError(f"seed must be an int, '{seed}' of type '{type(seed)}' was passed")
        super(RandomState, self).__init__(rd.MT19937(rd.SeedSequence(seed)))
