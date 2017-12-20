"""Hydrograph separation methods

Hydrograph separation refers to the process of decomposing streamflow series
into components representing water incoming from different sources. These are
generally assumed to be the quickflow and the baseflow.

The quickflow refers to the amount of water that reaches the stream quickly
by the means surface overland flow or near stream subsurface flow.

The baseflow refers to the amount of delayed water coming from groundwater
discharge.

"""

import numpy as np
import pandas as pd
from scipy.signal import lfilter, lfilter_zi, argrelmin, lfiltic


def eckhardt(q, alpha, beta, init_value=None, window=7):
    """Eckhardt recursive filter for hydrograph separation.

    Implementation of the Eckhardt digital recursive filter.

    .. math:
        q_{b(i)} = \frac{(1-\beta)\alpha q_{b(i-1)} + (1-\alpha)\beta q_i}
        {1-\alpha \beta}


    Parameters
    ----------
    q : array_like
        Time series of streamflow values.
    alpha : float
        Recession coefficient representing the proportion of remaining
        streamflow on the next time step. It should be strictly greather than
        0 and lower than 1.
    beta : float
        Maximum baseflow index representing the long term ratio between
        baseflow and streamflow.
    init_value : None, float, optional
        Initial value of the baseflow used to initiate Eckhardt filter. If None,
        then the first minimum of the first window
    window : int, optional


    Returns
    -------
    qb : array_like
        Time series of baseflow values.


    Examples
    --------

    >>> import numpy as np
    >>> q = np.array( [2.,2.1,2.4,2.8,2.6,2.4,2.3,2.2,2.1,2.,1.9])
    >>> baseflow = eckhardt(q, alpha=0.98, beta=0.66)
    >>> print(baseflow)

    References
    ----------

    Eckhardt, K. "How to Construct Recursive Digital Filters for Baseflow
    Separation." Hydrological Processes 19, no. 2 (February 15, 2005): 507-15.
    https://doi.org/10.1002/hyp.5675.

    """

    if (alpha <= 0) or (alpha >= 1.):
        raise ValueError(
            'Parameter alpha should be between range 0 < alpha < 1.')

    if (beta <= 0) or (beta >= 1.):
        raise ValueError(
            'Parameter beta should be between range 0 < beta < 1.')

    # Filter parameters are refactored to fit the scipy linear filter
    # parameters

    # input coefficient for the multiple backwards operators
    b = np.array([(1 - alpha) * beta / (1 - alpha * beta)])
    # output coefficient for the multiple backwards operators
    a = np.array([1, -(1 - beta) * alpha / (1 - alpha * beta)])

    # setting initial value

    if init_value is not None:
        qb0 = lfiltic(b, a, [init_value])
    else:
        locmin = min(q[:window])
        qb0 = lfiltic(b, a, [locmin])

    qb, _ = lfilter(b, a, q, zi=qb0)

    # find where qb is higher than q

    invalid = np.where(qb > q)

    # replace qb > q values by q

    qb[invalid] = q.iloc[invalid]

    return qb


def naive_eckhardt(q, alpha, beta):
    """naive implementation"""

    qb = np.zeros(len(q) + 1)
    qb[0] = q.iloc[:7].min()

    for i, qi in enumerate(q):
        factor1 = (1 - beta) * alpha * qb[i]
        factor2 = (1 - alpha) * beta * qi
        factor3 = 1 - alpha * beta
        qb[i + 1] = (factor1 + factor2) / factor3

    # strip initial value
    qb = qb[1:]

    # find where qb is higher than q

    invalid = np.where(qb > q)

    # replace qb > q values by q

    qb[invalid] = q.iloc[invalid]

    return qb


if __name__ == '__name__':
    import doctest

    doctest.testmod()
