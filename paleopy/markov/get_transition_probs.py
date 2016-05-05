
def get_transition_probs(x, classes=None):
    import numpy as np
    import pandas as pd
    """
    x is a dataframe or a time-series

    """
    ### ==============================================================================================================
    ### define list_ind as the list of unique ID if not defined
    if not classes:
        classes = np.unique(x.values)
    ### ==============================================================================================================
    transitions_matrix = []
    for cl in classes:
        list_mat = []
        list_freq = []
        for yi in np.unique(x.index.year):
            ### select the year
            mat = x.ix[str(yi),].values.flatten()
            ### select the class cl index + 1
            i = np.where(mat == cl)[0] + 1
            ### make sure we clip
            i = np.delete(i, np.where(i >= len(mat)))
            ### count the total number of days following the cluster "cl"
            list_freq.append(len(i))
            mat = mat[i,].tolist()
            ### count the number of occurences of each of the N clusters
            list_mat.append([mat.count(c) for c in classes])
        list_mat = np.array(list_mat)
        list_mat = list_mat.sum(0) * 1.0
        list_mat = list_mat / np.array(list_freq).sum()
        transitions_matrix.append(list_mat)
    transitions_matrix = np.array(transitions_matrix)
    return classes, transitions_matrix
