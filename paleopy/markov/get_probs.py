
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
