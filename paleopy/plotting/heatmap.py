def heatmap(df_anoms, df_probs, df_probs_MC): 
    sig90 = make_sig(df_probs, df_probs_MC, sig=90)
    sig95 = make_sig(df_probs, df_probs_MC, sig=95)
    sig99 = make_sig(df_probs, df_probs_MC, sig=99)
    
    df_anoms = df_anoms * 100
    df_anoms = df_anoms.T
    anoms = df_anoms.values[::-1,:]
    
    vmax = np.ceil(np.abs([anoms.min(),anoms.max()]).max())
    vmin = -1 * vmax

    width, height = plt.figaspect(anoms)
    fig, ax = plt.subplots(figsize=(10, 10))

    #fig = plt.figure(figsize=(12,10))
    #ax = fig.add_subplot(gs[2::,2::])

    X = np.arange(df_anoms.shape[1])
    Y = np.arange(0.5, df_anoms.shape[0] + 0.5)
    # im = ax.pcolormesh(X, Y, anoms,vmin=vmin, vmax=vmax, cmap=plt.cm.RdBu_r) 
    im = ax.pcolormesh(anoms,vmin=vmin, vmax=vmax, cmap=plt.cm.RdBu_r) 

    for i in range(anoms.shape[0]):
        for j in range(anoms.shape[1]): # maybe not hard code that 
            if sig95[i,j]:
                ax.plot(j+0.50,i+0.50, "s",  color='k', markersize=11, markerfacecolor='None')
            if sig99[i,j]:
                ax.plot(j+0.50,i+0.50, "+",  color='k', markersize=11, markerfacecolor='None')

    ax.set_yticks(np.arange(0.5, anoms.shape[0] + 0.5))
    ax.set_yticklabels(df_anoms.index[::-1])

    ax.set_xticks(np.arange(0.5, anoms.shape[1] + 0.5))
    ax.set_xticklabels(df_anoms.columns)


    cb = plt.colorbar(im, ticks=np.linspace(vmin, vmax, 11), boundaries=np.linspace(vmin, vmax, 11), drawedges=True)

    return fig