def plot_season_ts(p):
    r"""
    plots a seasonal time-series for a proxy
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from dateutil.relativedelta import relativedelta


    f, ax = plt.subplots(figsize=(8,5))

    if p.calc_anoms:

        y = p.ts_seas.loc[:,'anomalies']
        vmin = y.min() -0.1 * (np.abs(y.min()))
        vmax = y.max() +0.1 * (np.abs(y.min()))

        yd = p.ts_seas.loc[:,'d_anomalies']

        if p.detrend:
            ya = p.analogs.loc[:,'d_anomalies']
        else:
            ya = p.analogs.loc[:,'anomalies']

    else:

        y = p.ts_seas.loc[:,p.variable]
        vmin = y.min() -0.01 * (np.abs(y.min()))
        vmax = y.max() +0.01 * (np.abs(y.min()))

        yd = p.ts_seas.loc[:,'d_' + p.variable]

        if p.detrend:
            ya = p.analogs.loc[:,'d_' + p.variable]
        else:
            ya = p.analogs.loc[:,p.variable]

    ax.plot(y.index, y.values, 'steelblue', lw=2, label='{}'.format(p.variable))
    ax.plot(yd.index, yd.values, color='k', lw=2, label='{} (detrended)'.format(p.variable))
    ax.plot(ya.index, ya.values, 'ro', label='analog years')
    ax.vlines(ya.index, vmin, vmax, lw=5, alpha=0.3, label="")

    ax.set_xlim(y.index[0] - relativedelta(years=1), y.index[-1] + relativedelta(years=1))

    ax.set_ylim(vmin, vmax)
    ax.set_ylabel(p.variable +": " + p.dset_dict['units'], fontsize=14)

    ax.legend(framealpha=0.4, loc='best')

    if not p.qualitative and p.method == 'quintiles':
        [ax.axhline(b, color='magenta', zorder=1, alpha=0.5) for b in p.quintiles[1:-1]]

    # add a zero line if we deal with anomalies
    if p.calc_anoms:
        ax.axhline(0, color='k', linewidth=0.5)

    ax.grid()

    lyears = ",".join(map(str, p.analogs.index.year.tolist()))

    ax.set_title("Analog seasons for {} {} from {} {}:\n{}".format(p.season, p.sitename, p.dataset,
                                                                p.variable, lyears, fontsize=14))

    [l.set_fontsize(14) for l in ax.xaxis.get_ticklabels()]
    [l.set_fontsize(14) for l in ax.yaxis.get_ticklabels()]

    return f
