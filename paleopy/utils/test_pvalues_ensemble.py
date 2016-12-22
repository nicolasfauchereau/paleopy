def test_pvalues_ensemble(x, test=0.1, method='all'):
    """
    in the case of an ensemble of heteroclite proxies,
    test where ALL, ANY or a given percentage of the proxies
    composite anomalies are significant at the given signficance level

    needs to be apply as

    np.apply_along_axis(test_pvalues_ensemble, 0, xarray.data, **{'test':0.1, 'method':50})

    parameters
    ==========

    x = a one dimensional array

    """
    import numpy as np
    if method == 'all':
        ptest = np.all(x <= test).astype(int)
    elif method == 'any':
        ptest = np.any(x <= test).astype(int)
    else:
        try:
            percent = float(method) / 100.
            ptest = ((x.sum() / len(x)) >= percent).astype(int)
        except:
            print("cannot coerce method to be a float (percentage)")
    return ptest
