def get_probs(x, classes):
    import numpy as np
    """
    Returns the normalized frequencies (i.e. probability
            of occurrences) of each class (defined in classes)
            in x
    Parameters
    """
    if type(x) != list:
        try:
            x = x.flatten().tolist()
        except:
            x = x.values.flatten().tolist()
    if type(classes) != list:
        classes = classes.tolist()
    return np.array([float(x.count(c)) for c in classes]) / len(x)

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
            ### count the number of occurences of each of the 12 clusters
            ### TODO (function with N = number of clusters)
            list_mat.append([mat.count(c) for c in classes])
        list_mat = np.array(list_mat)
        list_mat = list_mat.sum(0) * 1.0
        list_mat = list_mat / np.array(list_freq).sum()
        transitions_matrix.append(list_mat)
    transitions_matrix = np.array(transitions_matrix)
    return classes, transitions_matrix

def simulate(mu,p,N,states):
    """
    simulate a discrete K states markov chain time-series given:
        INPUTS
    1) mu (np.array) = the probability of each discrete state, must be of length K (number of states)
    2) p (np.array) = the transition matrix, describing the probability of transitions
    from each state to every other
    must be of shape (K,K) with first row describing probablity of transition
    from state 0 to states 0(i.e. itself) to K (number of states)
    3) N (integer) = the lenth of the one dimensional time-series to simulate
    4) states = a list of states (can ve any type)
    returns:
        OUTPUTS
    x = a np.array of length N with simulated states (from 0 to K-1)
    """
    states = dict(enumerate(states))

    x = np.zeros((N,)).astype(np.int32)

    def rMC(mu):
        """
        Draw an integer in 0 -> len(mu)
        given the probabilities in mu
        """
        import numpy as np
        u = np.random.random_sample()
        i = 0.
        s = mu[0]

        while ( (u > s) & (i < len(mu))  ):
            i += 1
            s = s + mu[i]
        index = i
        return index

    """
    Initialize the sequence
    """
    x[0] = rMC(mu)

    """
    Now iterate, at each time step, the class (i.e. weather regime)
    is drawn according to the probability of transitions defined
    in p
    """
    for i in range(N-1):
        x[i+1] = rMC(p[x[i]])

    x = x.astype(np.int32)

    """
    Get actual values (i.e. strings) from keys in dictionnary
    """
    return np.asarray([states[k] for k in x])

def simulate_2D(states, freqs, trans, N=1000, P=1000):
	"""

	"""
    import numpy as np
	x = np.empty((N,P),dtype=np.int32)

	x[0,:] = np.random.choice(len(states), size=(P,), replace=True, p=freqs)

	for j in range(P):
		for i in range(N-1):
			x[i+1,j] = np.random.choice(len(states), size=1, replace=True, p=trans[x[i,j],:])

	return dict(enumerate(states)), x
