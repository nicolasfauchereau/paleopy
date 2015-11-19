def make_sig(df_probs, df_probs_MC, sig=90, array=True):
    import pandas as pd
    sigdf = []
    for icol in range(df_probs.shape[1]):
        col = df_probs.iloc[:,icol]
        sigdf.append((col < df_probs_MC.loc[:,str((100 - sig))]) | (col > df_probs_MC.loc[:,str(sig)]))
    sigdf = pd.concat(sigdf, axis=1)
    sigdf.columns = df_probs.columns
    if array:
        sigdf = sigdf.T
        sigarray = sigdf.values[::-1,:]
        return sigarray
    else:
        return sigdf
