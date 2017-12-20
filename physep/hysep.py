# coding=utf-8
import numpy as np
import pandas as pd
from scipy.ndimage.filters import minimum_filter1d, generic_filter
from scipy.ndimage.measurements import label
from scipy.signal import argrelextrema


def minimum_filter(ts, **kwargs):
    """Return stationary base flow
    
    The base flow is set to the minimum observed flow.
    
    :param ts: 
    :return: 
    """
    minimum = min(ts)
    out_values = minimum * np.ones(len(ts))
    baseflow = pd.Series(data=out_values, index=ts.index)
    quickflow = ts - baseflow
    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'
    return baseflow, quickflow


def fixed_interval_filter(ts, size):
    """USGS HYSEP fixed interval method
    
    The USGS HYSEP fixed interval method as described in `Sloto & Crouse, 1996`_.
    
    .. _Slot & Crouse, 1996:
        Sloto, Ronald A., and Michele Y. Crouse. “HYSEP: A Computer Program for Streamflow Hydrograph Separation and 
        Analysis.” USGS Numbered Series. Water-Resources Investigations Report. Geological Survey (U.S.), 1996. 
        http://pubs.er.usgs.gov/publication/wri964040.

    
    :param size: 
    :param ts: 
    :return: 
    """
    intervals = np.arange(len(ts)) // size
    baseflow = pd.Series(data=ts.groupby(intervals).transform('min'), index=ts.index)
    quickflow = ts - baseflow

    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow


def sliding_interval_filter(ts, size):
    """USGS HYSEP sliding interval method
    
        The USGS HYSEP sliding interval method as described in `Sloto & Crouse, 1996`_.
        
        The flow series is filter with scipy.ndimage.genericfilter1D using numpy.nanmin function
        over a window of size `size`
    
    .. _Slot & Crouse, 1996:
        Sloto, Ronald A., and Michele Y. Crouse. “HYSEP: A Computer Program for Streamflow Hydrograph Separation and 
        Analysis.” USGS Numbered Series. Water-Resources Investigations Report. Geological Survey (U.S.), 1996. 
        http://pubs.er.usgs.gov/publication/wri964040.
    
    :param size: 
    :param ts: 
    :return: 
    """
    # TODO ckeck the presence of nodata
    if (ts.isnull()).any():
        blocks, nfeatures = label(~ts.isnull())
        block_list = [ts[blocks == i] for i in range(1, nfeatures + 1)]
        na_df = ts[blocks == 0]
        block_bf = [pd.Series(data=minimum_filter1d(block, size, mode='reflect'), index=block.index) for block in
                    block_list]
        baseflow = pd.concat(block_bf + [na_df], axis=0)
        baseflow.sort_index(inplace=True)
    else:
        baseflow = pd.Series(data=minimum_filter1d(ts, size, mode='reflect'), index=ts.index)

    quickflow = ts - baseflow

    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow


def local_minimum_filter(ts, size):
    """USGS HYSEP local minimum method
    
        The USGS HYSEP local minimum method as described in `Sloto & Crouse, 1996`_.
    
    .. _Slot & Crouse, 1996:
        Sloto, Ronald A., and Michele Y. Crouse. “HYSEP: A Computer Program for Streamflow Hydrograph Separation and 
        Analysis.” USGS Numbered Series. Water-Resources Investigations Report. Geological Survey (U.S.), 1996. 
        http://pubs.er.usgs.gov/publication/wri964040.
    
    :param size: 
    :param ts: 
    :return: 
    """

    origin = int(size) / 2
    baseflow_min = pd.Series(generic_filter(ts, _local_minimum, footprint=np.ones(size)), index=ts.index)
    baseflow = baseflow_min.interpolate(method='linear')
    # interpolation between values may lead to baseflow > streamflow
    errors = (baseflow > ts)
    while errors.any():
        print 'hello world'
        error_labelled, n_features = label(errors)
        error_blocks = [ts[error_labelled == i] for i in range(1, n_features + 1)]
        error_local_min = [argrelextrema(e.values, np.less)[0] for e in error_blocks]
        print error_local_min
        break
    quickflow = ts - baseflow
    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow


def eckhardt_filter(ts, size):
    return fixed_interval_filter(ts, size)

def _local_minimum(window):
    win_center_ix = len(window) / 2
    win_center_val = window[win_center_ix]
    win_minimum = np.min(window)
    if win_center_val == win_minimum:
        return win_center_val
    else:
        return np.nan


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    data = pd.Series.from_csv('test/hysep_data.csv',
                              parse_dates=True,
                              index_col=0,
                              header=0)
    b, q = local_minimum_filter(data, 93)
    plt.plot(data)
    plt.plot(b)
    plt.show()
