"""
A collection of simulation functions for State Space Models.
"""
import numpy as np
from .random import RandomState

__all__ = ["generate_ssm", "generate_pendulum"]


def _atleast2d(*args):
    return tuple(np.atleast_2d(elem) for elem in args)


def generate_ssm(m_0, A, Q, H, R, steps, random_state):
    """Samples from a state space model given parameters and a random state

    Parameters
    ----------
    m_0 : (M,) array_like
        Initial mean of the state
    P_0 : (M, M) array_like
        Initial covariance of the state
    A : (M, M) or (M, M) array_like
        Transition matrix
    Q : (M, M) array_like
        Transition covariance
    H : (M, N) array_like
        Observation matrix
    R : (N, N) array_like
        Observation covariance
    steps : int
        Number of steps simulated
    random_state : RandomState
        Random state used for pseudo-random numbers generation

    Returns
    -------
    states : (steps, M) ndarray
        The true states
    observations : (steps, N) ndarray
        The noisy observations

    Examples
    --------
    >>> M, N = 2, 1
    >>> m_0 = np.zeros(M)
    >>> P_0 = Q  = [[0.4, -0.2],
    ...             [-0.2, 0.5]]
    >>> A = np.zeros((M, M))
    >>> H = [0., 0.]
    >>> R = 0.5
    >>> states, observations = generate_ssm(m_0, P_0, A, Q, H, R, 10000, RandomState(5))
    >>> est_cov = np.cov(states, rowvar=False)
    >>> est_error = np.cov(observations, rowvar=False)
    >>> cov_close = np.allclose(est_cov, Q, atol=1e-2)
    >>> error_close = np.allclose(est_error, R, atol=1e-2)
    >>> cov_close & error_close
    True
    """
    if not isinstance(random_state, RandomState):
        raise TypeError(f"random_state must be an instance of {RandomState}, "
                        f"'{random_state}' of type '{type(random_state)}' was given")

    m_0 = np.atleast_1d(m_0)
    A, Q, H, R = _atleast2d(A, Q, H, R)

    M = m_0.shape[-1]
    N = R.shape[-1]
    states = np.empty((steps, M))
    observations = np.empty((steps, N))

    chol_Q = np.linalg.cholesky(Q)
    chol_R = np.linalg.cholesky(R)

    state = m_0
    for i in range(steps):
        state = A @ state + chol_Q @ random_state.randn(M)
        states[i, :] = state
        obs = H @ state + chol_R @ random_state.randn(N)
        observations[i, :] = obs

    return states, observations


def generate_pendulum(m_0, g, Q, dt, R, steps, random_state, cluttered_probability=0., clutter_range=(-2, 2)):
    """ Samples from a noisy pendulum submitted to a gravitational pull g a random state.
    The state represents the angle and the angular moment, the measurement is the sine of the angle:
    the horizontal position of the pendulum.

    Parameters
    ----------
    m_0 : (2,) array_like
        Initial mean of the state
    g : float
        Gravitational pull (g-force) in N/kg (earth is ~9.81)
    Q : (2, 2) array_like
        Transition covariance coming from the discretisation of the model
    dt : float
        Time between each measurement
    R : float
        Observation variance
    steps : int
        Number of steps simulated
    random_state : RandomState
        Random state used for pseudo-random numbers generation
    cluttered_probability : float, optional
        What are the chances that the observations are cluttered
    clutter_range: tuple of float
        When observation is cluttered, it's replaced by a uniform in this range


    Returns
    -------
    timeline: (steps) ndarray
        The observation times
    states : (steps, M) ndarray
        The true states
    observations : (steps, N) ndarray
        The noisy observations
    """
    if not isinstance(random_state, RandomState):
        raise TypeError(f"random_state must be an instance of {RandomState}, "
                        f"'{random_state}' of type '{type(random_state)}' was given")

    m_0 = np.atleast_1d(m_0)
    Q = np.atleast_2d(Q)

    states = np.empty((steps, 2))
    observations = np.empty(steps)

    chol_Q = np.linalg.cholesky(Q)
    sqrt_R = np.sqrt(R)

    state = m_0

    for i in range(steps):
        state = np.array([state[0] + dt * state[1],
                          state[1] - g * dt * np.sin(state[0])])
        state = state + chol_Q @ random_state.randn(2)
        states[i, :] = state

        observations[i] = np.sin(state[0]) + sqrt_R * random_state.randn()

    if cluttered_probability > 0:
        cluttered_ind = random_state.rand(steps) < cluttered_probability
        clutter_multiplier = clutter_range[1] - clutter_range[0]
        observations[cluttered_ind] = random_state.rand(
            cluttered_ind.astype(np.int_).sum()) * clutter_multiplier - clutter_multiplier / 2.

    return np.arange(dt, (steps + 1) * dt, dt), states, observations
