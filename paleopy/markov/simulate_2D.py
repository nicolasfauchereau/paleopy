import numpy as np

def simulate_2D(states, freqs, trans, N=1000, P=1000):
	"""
	simulates a P vectors of size N given a unique vector of states,
	a vector of frequencies (must sum to 1) and a transition matrix
	"""
	x = np.empty((N,P),dtype=np.int32)

	x[0,:] = np.random.choice(len(states), size=(P,), replace=True, p=freqs)

	for j in range(P):
		for i in range(N-1):
			x[i+1,j] = np.random.choice(len(states), size=1, replace=True, p=trans[x[i,j],:])

	return dict(enumerate(states)), x
