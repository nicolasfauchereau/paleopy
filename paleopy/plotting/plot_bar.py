def plot_bar(wr, ensemble=True):
    """
    """
    
    # if the climatological probabilities have 
    # not been calculated, call the `_get_clim_probs` 
    # method
    if not(hasattr(wr, 'clim_probs')): 
        clim_probs = wr._get_clim_probs()
    
    # if the anomalies and the observed 
    # as well as the MC probabilities have 
    # not been calculated, call the `probs_anomalies`
    # method
    if not(hasattr(wr, 'df_anoms')): 
        wr.probs_anomalies(kind='many')
        
    if wr.df_anoms.shape[1] < 5:
        print("""
        not enough proxies in the ensemble, size was {}
        need at list 5
        """.format(wr.df_anoms.shape[1]))
        raise Exception("size error")
            
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)
    ax1.bar(np.arange(0.5, len(clim_probs)+0.5), clim_probs * 100, color="0.8", width=1., alpha=0.8)
    ax1.set_ylabel("climatological frequency %")
    ax2 = ax1.twinx()
    ax2.axhline(0,color='k', lw=1, linestyle='dashed')
    ax2.set_ylabel("% change in frequency")
    
    df_anoms = wr.df_anoms * 100
    
    bp = df_anoms.T.boxplot(ax=ax2, patch_artist=True);
    
    wrone = copy(wr)
    wrone.probs_anomalies(kind='one')
    testb = (wrone.df_probs['ensemble'] < wrone.df_probs["10"]) \
    | (wrone.df_probs['ensemble'] > wrone.df_probs["90"]).values
    
    
    for i,b in enumerate(bp['boxes']):
        if wrone.df_anoms.iloc[i,:].values >= 0: 
            plt.setp(b,facecolor='r', edgecolor='r', alpha=0.5)
            plt.plot(i+1,wrone.df_anoms.iloc[i,:]*100, 'r*')
            if testb[i]: 
                plt.setp(b,facecolor='r', edgecolor='r', alpha=0.9)
        else: 
            plt.setp(b,facecolor='b', edgecolor='b', alpha=0.5)
            plt.plot(i+1,wrone.df_anoms.iloc[i,:]*100, 'b*')
            if testb[i]: 
                plt.setp(b,facecolor='b', edgecolor='r', alpha=0.9)
    
    plt.setp(bp['fliers'], color='gray', marker='+',visible=True)
    plt.setp(bp['medians'],color='k',linewidth=1)
    
    ax2.set_title(w.classification, fontsize=14)
    
    return fig