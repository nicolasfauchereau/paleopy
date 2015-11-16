def plot_bar(wr, ensemble=True)
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
            
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.bar(np.arange(0.5, len(clim_probs)+0.5), clim_probs * 100, color="0.8", width=1., alpha=0.8)
    ax1.set_ylabel("frequency %")
    ax2 = ax1.twinx()
    ax2.axhline(0,color='k', lw=1, linestyle='dashed')
    ax2.set_ylabel("% change in frequency")
    
    df_anoms = wr.df_anoms * 100
    
    bp = df_anoms.T.boxplot(ax=ax2, patch_artist=True);
    
    wrone = copy.copy(wr)
    wrone.probs_anomalies(kind='one')
    
    for i,b in enumerate(bp['boxes']):
        if wrone.df_anoms.iloc[i,:].values >= 0: 
            plt.setp(b,facecolor='red', alpha=0.5)
        else: 
            plt.setp(b,facecolor='blue', alpha=0.5)
    
    return wrone