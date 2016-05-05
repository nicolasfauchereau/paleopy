def simulate(mu,p,N,states):
    import numpy as np
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
