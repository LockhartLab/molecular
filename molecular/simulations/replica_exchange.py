
import numpy as np


def temp_schedule(temp_min=300, temp_max=440, n_temps=40, mode='geometric'):
    r"""
    Create a temperature schedule that could be used, for instance, with replica exchange.

    There are several choices for `mode`:
        * "geometric"

        .. math :: T_i = T_1 \frac{T_R}{T_1}^{\frac{i-1}{R-1}}

        * "linear"

        .. math :: T_i = T_1 + (i-1) \frac{T_R-T_1}{R-1}

        * "parabolic" (Note if `n_replicas` is even, `temp_max` won't directly be sampled).

        .. math :: T_i = T_1 - \frac{T_R-T_1}{\left( \frac{R-1}{2} \right) ^2} (i-1) (i-R)

    Parameters
    ----------
    temp_min : float
        Lowest temperature
    temp_max : float
        Highest temperature
    n_temps : int
        Number of temperatures
    mode : str
        Mode to produce schedule. Valid options include "geometric", "linear", "parabolic". (Default: "geometric")

    Returns
    -------
    numpy.ndarray
        Temperature schedule
    """

    mode = mode.lower()

    if mode in 'geometric':
        schedule = temp_min * np.power(temp_max / temp_min, np.arange(n_temps) / (n_temps - 1.), dtype=np.float64)

    elif mode in 'linear':
        schedule = np.linspace(start=temp_min, stop=temp_max, num=n_temps)

    elif mode in 'parabolic':
        temp_range = temp_max - temp_min
        temp_ind = np.arange(n_temps)
        schedule = temp_min - (temp_range / np.square((n_temps - 1.) / 2.)) * temp_ind * (temp_ind - n_temps + 1.)

    else:
        raise AttributeError(f'mode {mode} not supported')

    return schedule


if __name__ == '__main__':
    # print(reptemp(temp_min=300, temp_max=440, n_replicas=40, mode='geometric'))
    print(temp_schedule(temp_min=310, temp_max=500, n_temps=5, mode='parabolic'))