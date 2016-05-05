def bar(wr, sig):
    """
    """

    from copy import copy
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # if the climatological probabilities have
    # not been calculated, call the `_get_clim_probs`
    # method
    if not(hasattr(wr, 'clim_probs')):
        clim_probs = wr._get_clim_probs()

    if wr.parent.description == 'ensemble':

        # if the anomalies and the observed
        # as well as the MC probabilities have
        # not been calculated, call the `probs_anomalies`
        # method
        if not(hasattr(wr, 'df_anoms')):
            wr.probs_anomalies(kind='many')

        if wr.df_anoms.shape[1] < 3:
            print("""
            not enough proxies in the ensemble, size was {}
            need at list 5
            """.format(wr.df_anoms.shape[1]))
            raise Exception("size error")

        fig = plt.figure(figsize=(6,6))

        """
        AXES number 1: climatological frequencies
        """

        ax1 = fig.add_subplot(111)

        ax1.bar(np.arange(0.5, len(clim_probs)+0.5), clim_probs * 100, color="0.8", width=1., alpha=0.8)
        ax1.set_ylabel("climatological frequency %", fontsize=14)

        [l.set_fontsize(14) for l in ax1.yaxis.get_ticklabels()]
        [l.set_fontsize(14) for l in ax1.xaxis.get_ticklabels()]
        [l.set_rotation(90) for l in ax1.xaxis.get_ticklabels()]

        """
        AXES number 2: frequency anomalies
        """
        ax2 = ax1.twinx()
        ax2.axhline(0,color='k', lw=1, linestyle='dashed')
        ax2.set_ylabel("% change in frequency", fontsize=14)

        df_anoms = wr.df_anoms * 100

        bp = df_anoms.T.boxplot(ax=ax2, patch_artist=True);

        wrone = copy(wr)

        wrone.probs_anomalies(kind='one')

        testb = (wrone.df_probs['ensemble'] < wrone.df_probs["{}".format(sig)]) \
        | (wrone.df_probs['ensemble'] > wrone.df_probs["{}".format((100-sig))]).values

        for i,b in enumerate(bp['boxes']):
            if wrone.df_anoms.iloc[i,:].values >= 0:
                plt.setp(b,facecolor='coral', edgecolor='coral', alpha=0.5)
                plt.plot(i+1,wrone.df_anoms.iloc[i,:]*100, 'r*')
                if testb[i]:
                    plt.setp(b,facecolor='r', edgecolor='r', alpha=0.9)
            else:
                plt.setp(b,facecolor='steelblue', edgecolor='steelblue', alpha=0.5)
                plt.plot(i+1,wrone.df_anoms.iloc[i,:]*100, 'b*')
                if testb[i]:
                    plt.setp(b,facecolor='b', edgecolor='b', alpha=0.9)

        plt.setp(bp['fliers'], color='gray', marker='+',visible=True)
        plt.setp(bp['medians'],color='k',linewidth=1)

        [l.set_fontsize(14) for l in ax2.xaxis.get_ticklabels()]
        [l.set_fontsize(14) for l in ax2.yaxis.get_ticklabels()]
        [l.set_rotation(90) for l in ax2.xaxis.get_ticklabels()]
        ax2.set_title("{} Weather Regimes".format(wr.classification), fontsize=14)


    # if the object passed is NOT an ensemble, but a proxy object
    else:
        # get the anomalies
        if not(hasattr(wr, 'df_anoms')):
            wr.probs_anomalies(kind='one')

        df_anoms = wr.df_anoms * 100

        nregimes = len(df_anoms)

        testb = (wr.df_probs[wr.parent.sitename] < wr.df_probs["{}".format(sig)]) \
            | (wr.df_probs[wr.parent.sitename] > wr.df_probs["{}".format((100-sig))]).values

        # first subplot, climatological frequencies

        fig, axes = plt.subplots(nrows=2,ncols=1, figsize=(nregimes*0.8,9), sharex=True)
        axes = axes.flatten()
        fig.subplots_adjust(hspace=0)

        """
        first axes: climatological frequencies
        """
        ax1 = axes[0]

        ax1.bar(np.arange(0.5, nregimes+0.5), clim_probs * 100, color="0.8", width=1., alpha=0.8)

        ax1.set_ylabel("climatological frequency (%)", fontsize=14)

        [l.set_fontsize(13) for l in ax1.yaxis.get_ticklabels()]

        ax1.set_title("{} Weather Regimes".format(wr.classification), fontsize=14)

        ax1.yaxis.grid(True)
        ax1.set_axisbelow(True)

        [ax1.axvline(l, color='k', linestyle=':', zorder=-1) for l in np.arange(1.5, nregimes)]


        """
        second axes: anomalies
        """
        ax2 = axes[1]

        # loop over the dataframe containing the anomalies, checks the
        # test, and plot ...
        for i in range(nregimes):
            if df_anoms.iloc[i,:].values >= 0:
                if testb.iloc[i]:
                    ax2.bar(i+0.5, df_anoms.iloc[i,:].values, color='r',width=1., alpha=1.)
                else:
                    ax2.bar(i+0.5, df_anoms.iloc[i,:].values, color='coral',width=1., alpha=.3)
            else:
                if testb.iloc[i]:
                    ax2.bar(i+0.5, df_anoms.iloc[i,:].values, color='b',width=1., alpha=1.)
                else:
                    ax2.bar(i+0.5, df_anoms.iloc[i,:].values, color='steelblue',width=1., alpha=.3)

        # move the ticks and labels to the right
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")

        ax2.set_xticks(range(1, nregimes+1))
        ax2.set_xticklabels(df_anoms.index.tolist())
        ax2.set_xlim(0.5, nregimes + 0.5)

        ax2.set_ylabel("change in frequency (%)", fontsize=14)

        ax2.set_xlabel("weather regime", fontsize=14)

        [l.set_fontsize(13) for l in ax2.yaxis.get_ticklabels()]

        [l.set_fontsize(14) for l in ax2.xaxis.get_ticklabels()]

        ax2.yaxis.grid(True)
        ax2.set_axisbelow(True)

        [ax2.axvline(l, color='k', linestyle=':', zorder=-1) for l in np.arange(1.5, nregimes)]

    return fig
