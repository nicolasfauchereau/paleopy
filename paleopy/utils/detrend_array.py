def detrend_array(x):
    import numpy as np
    from matplotlib.mlab import detrend_linear
    """
    detrend function to apply to a xarray.DataArray

    usage:

    dset['detrended_var'] = (('time', 'latitudes', 'longitudes'), np.apply_along_axis(detrend_array, 0, dset['var'].data))
    """
    if any(np.isnan(x)):
        z = np.ones(len(x)) * np.nan
    else:
        z = detrend_linear(x) + x.mean()
    return z
