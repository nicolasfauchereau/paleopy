def heatmap(wr):
    """
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from ..utils import make_sig

    df_anoms = wr.df_anoms
    df_probs = wr.df_probs
    df_probs_MC = wr.df_probs_MC
    classification = wr.classification

    sig90 = make_sig(df_probs, df_probs_MC, sig=90)
    sig95 = make_sig(df_probs, df_probs_MC, sig=95)
    sig99 = make_sig(df_probs, df_probs_MC, sig=99)

    df_anoms = df_anoms * 100
    df_anoms = df_anoms.T
    anoms = df_anoms.values[::-1,:]

    vmax = np.ceil(np.abs([anoms.min(),anoms.max()]).max())
    vmin = -1 * vmax

    width, height = plt.figaspect(anoms)

    fig, ax = plt.subplots(figsize=(10,10))
    fig.subplots_adjust(left=0.2, right=0.99)


    X = np.arange(df_anoms.shape[1])
    Y = np.arange(0.5, df_anoms.shape[0] + 0.5)

    # im = ax.pcolormesh(X, Y, anoms,vmin=vmin, vmax=vmax, cmap=plt.cm.RdBu_r)
    im = ax.pcolormesh(anoms,vmin=vmin, vmax=vmax, cmap=plt.cm.RdBu_r)

    # this is a small adjustment to align the vertical lines to the
    # cells boundaries
#     adj = 0.02
    adj = 0

    for i in range(anoms.shape[0]):
        for j in range(anoms.shape[1]): # maybe not hard code that
            if sig95[i,j]:
                ax.plot(j+0.50,i+0.50, "s",  color='k', markersize=11, \
                        markerfacecolor='None', zorder=1, mew=1.25)
            if sig99[i,j]:
                ax.plot(j+0.50+adj,i+0.50-adj, "x",  color='k', markersize=10, \
                        markerfacecolor='None', zorder=1, mew=1.5)

    ax.set_yticks(np.arange(0.5, anoms.shape[0] + 0.5))
    ax.set_yticklabels(df_anoms.index[::-1])

    ax.set_xticks(np.arange(0.5, anoms.shape[1] + 0.5))
    ax.set_xticklabels(df_anoms.columns)

    [l.set_fontsize(12) for l in ax.yaxis.get_ticklabels()]
    [l.set_fontsize(12) for l in ax.xaxis.get_ticklabels()]
    [l.set_rotation(90) for l in ax.xaxis.get_ticklabels()]

    # this is a small adjustment to align the vertical lines to the
    # cells boundaries
    #     adj = 0.01
    adj = 0

    ax.hlines(np.arange(1,anoms.shape[0]), 0, anoms.shape[1])
    ax.vlines(np.arange(1-adj,anoms.shape[1]-adj), 0, anoms.shape[0], color='0.1')

    # remobes the ticks on each axis
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # plots the colorbar
    cb = plt.colorbar(im, ticks=np.linspace(vmin, vmax, 11), boundaries=np.linspace(vmin, vmax, 11), drawedges=True)
    cb.set_label('% change in frequency', fontsize=14)
    [l.set_fontsize(12) for l in cb.ax.yaxis.get_ticklabels()]

    ax.set_title("{} Weather Regimes".format(classification), fontsize=14)

    return fig
