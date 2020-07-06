"""
A stats utility file used to compute usual statistics
"""
import numpy as np

__all__ = ["rmse"]


def rmse(x, y):
    """Computes the RMSE between x and y along all their dimensions, x and y must have the same dimensions.

        Parameters
        ----------
        x : (L,...) array_like
            Initial mean of the state
        y : (L, ...) array_like
            Initial covariance of the state

        Returns
        -------
        out : float
            The RMSE
    """
    return np.sqrt(np.mean(np.sum(np.square(x - y), -1)))
