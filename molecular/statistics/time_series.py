"""
time_series.py

author: C. Lockhart <chris@lockhartlab.org>
language: Python3
"""

from numba import njit
import numpy as np


# noinspection PyShadowingNames
def acorr(a):
    r"""
    Compute the lagged autocorrelation of an observation :math:`a`.

    .. math :: \rho(\tau) = \frac{\gamma(\tau)}{\gamma(0)}

    The function :math:`\gamma` is from :ref:`acov`.

    Parameters
    ----------
    a : numpy.ndarray
        1D array.

    Returns
    -------
    numpy.ndarray
        Autocorrelation function
    """

    # Compute autocovariance
    gamma = acov(a)
    return gamma / gamma[0]


# noinspection PyShadowingNames
def _acorr(a):
    gamma = _acov(a)
    return gamma / gamma[0]


def _acorr_test(a, decimal=7, plot=True):
    rho0 = acorr(a)
    rho1 = _acorr(a)
    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        x = np.arange(len(rho0))
        plt.plot(x, rho0)
        plt.plot(x, rho1)
        plt.show()
    np.testing.assert_almost_equal(rho0, rho1, decimal=decimal)


# noinspection DuplicatedCode,PyShadowingNames
def acov(a, fft=False):
    r"""
    Compute the unbiased auto-covariance function :math:`\gamma(k)` from dataset :math:`a` with :math:`N` samples for
    lag-time :math:`k`.

    .. math :: \gamma(k)=\frac{1}{N-k}\sum_{t=1}^{N-k}(a_t - \mu)(a_{t+k} - \mu)

    Here, :math:`\mu=\frac{1}{N}\sum_t^Na_t`.

    Note: by default, all :math:`k` from 0 to N-1 are evaluated. Sampling deteriorates rapidly as :math:`k` increases.
    There is also a *biased* estimator for the autocovariance, which changes the denominator from :math:`n-k` to
    :math:`n` and has an effect of reducing the fluctuations due to error at large :math:`k`. To compute this, see
    :ref:`statsmodels.tsa.stattools.acf`.

    Parameters
    ----------
    a : numpy.ndarray
        1D array.
    fft : bool
        Should Fast-Fourier Transform be used?

    Returns
    -------
    numpy.ndarray
        Autocovariance function
    """

    # We only want to compute this for 1D arrays
    if a.ndim != 1:
        raise AttributeError('must be 1D')

    # Remove the mean from the observations
    a -= np.mean(a)

    # Compute autocovariance gamma
    if not fft:
        # Compute a * a for all lag times using correlate and the unbiased autocovariance
        len_a = len(a)
        gamma = np.correlate(a, a, mode='full')[(len_a - 1):] / (len_a - np.arange(len_a))

    else:
        # Use statsmodels to compute FFT
        # noinspection PyPackageRequirements
        from statsmodels.tsa.stattools import acovf
        gamma = acovf(a, adjusted=True, demean=False, fft=True, missing='none', nlag=None)

    # Return
    return gamma


@njit
def _acov(a):
    n = len(a)
    u = np.mean(a)
    gamma = np.zeros(n)
    for k in range(n):
        gamma[k] = np.mean((a[:n - k] - u) * (a[k:] - u))
    return gamma


def _acov_test(a):
    import time
    start_time = time.time()
    gamma0 = acov(a)
    end_time = time.time()
    print('acov={}'.format(end_time - start_time))

    start_time = time.time()
    gamma1 = acov(a, fft=True)
    end_time = time.time()
    print('acov={}'.format(end_time - start_time))

    start_time = time.time()
    gamma2 = _acov(a)
    end_time = time.time()
    print('acov={}'.format(end_time - start_time))

    np.testing.assert_almost_equal(gamma0, gamma1)
    np.testing.assert_almost_equal(gamma0, gamma2)


# noinspection PyShadowingNames
def inefficiency(a):
    """
    Compute the statistical inefficiency :math:`g` from the equilibration time \tau_{eq}.

    .. math :: g = 1 + 2\tau_{eq}.

    Parameters
    ----------
    a : numpy.ndarray

    Returns
    -------
    float
        Statistical inefficiency.
    """

    return 1. + 2. * teq(a)


# Estimate the standard error from the correlation time
# noinspection PyShadowingNames
def sem_tcorr(a, tol=1e-3):
    """
    Estimate of standard error of the mean derived from the correlation time.

    The main assumption is that sqrt(N) ~ sqrt(len(a) / tcorr(a))

    Should only be used for continuous time series data, e.g., from molecular dynamics. Discontinuous trajectories as
    produced by replica exchange or Monte Carlo are not applicable.

    Parameters
    ----------
    a : numpy.ndarray
    tol : float

    Returns
    -------

    """

    return np.std(a) * np.sqrt(tcorr(a, tol=tol) / len(a))


# TODO how is this different than teq?
# noinspection PyShadowingNames
def tcorr(a, tol=1e-3):
    """
    Compute correlation time from the autocorrelation function.

    Parameters
    ----------
    a : numpy.ndarray
    tol : float

    Returns
    -------
    int
        Index of `a` where first converged, i.e., always less than `tol`.
    """

    # Autocorrelation function
    # noinspection PyShadowingNames
    rho = acorr(a)

    # Accumulate the maximum for all time points
    # Strictly, this is performed on the reversed and absolute autocorrelation function rho (i.e., looking backwards)
    rho_rev_max = np.maximum.accumulate(np.abs(rho[::-1]))

    # Find rho_rev_max < tol
    rho_rev_equil = rho_rev_max < tol

    # If there are no instances where we are equilibrated,
    if np.sum(rho_rev_equil) == 0:
        raise ValueError(f'autocorrelation not converged with tol={tol}')

    # Return N - first observed convergence
    return len(rho) - np.argmin(rho_rev_equil)


# noinspection PyShadowingNames
def _tcorr(a, tol=1e-3):
    """
    Compute correlation time from the autocorrelation function.

    Parameters
    ----------
    a
    tol

    Returns
    -------

    """

    rho = acorr(a)

    for t in np.arange(len(rho)):
        if np.max(np.abs(rho[t:])) < tol:
            return t

    raise AttributeError(f'no correlation time found with tol={tol}')


def _tcorr_test(a, tol=1e-3):
    tau0 = tcorr(a, tol=tol)
    tau1 = _tcorr(a, tol=tol)
    np.testing.assert_equal(tau0, tau1)


# noinspection PyShadowingNames
def teq(a):
    r"""
    Compute equilibration time :math:`\tau_{eq}` based on the autocorrelation function :math:`\rho(t)`.

    .. math :: \tau_{eq} = \sum_{t=1}^{T} (1 - \frac{t}{T}) \rho_{t}

    Parameters
    ----------
    a : numpy.ndarray

    Returns
    -------
    float
        Equilibration time.
    """

    rho_ = acorr(a)[1:]
    rho = rho_[:np.min(np.where(rho_ < 0))]

    t_max = len(a)
    t = np.arange(1, len(rho) + 1)

    return np.sum((1. - t / t_max) * rho)


if __name__ == '__main__':
    n = 100000
    a = np.random.normal(loc=0, scale=10, size=n)
    print(teq(a))
    # _acov_test(a)
    # _acorr_test(a, decimal=5)
    # _tcorr_test(a)

    # print(np.std(a))
    # print(np.sqrt(np.cov(a)))
    # print(np.sqrt(xx - x*x))
    # print(correlation_time(a))
    # print(sem_tcorr(a))

    # import time
    # start_time = time.time()
    # rho = acorr(a)
    # end_time = time.time()
    # print(end_time - start_time)
    # start_time = time.time()
    # tau = tcorr(a)
    # end_time = time.time()
    # print(end_time - start_time)
