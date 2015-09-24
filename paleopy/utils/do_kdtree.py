def do_kdtree(lons, lats, point):
    """
    find the nearest point in [lons, lats] to 'point'
    """
    import numpy as np
    from scipy.spatial import cKDTree
    xy = np.dstack((lons, lats))[0]
    dtree = cKDTree(xy)
    dist, ind = dtree.query(point)
    return xy[ind]
