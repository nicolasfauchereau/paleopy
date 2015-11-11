def haversine(coord1, coord2):
    from math import radians, cos, sin, asin, sqrt
    """ return km
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    Parameters
    ----------
        coord1: coordinates (lon, lat) of the first location,
        coord1: coordinates (lon, lat) of the second location,

    Returns
    -------
        km: the distance between location 1 and 2 in km
    """
    # convert decimal degrees to radians
    lon1, lat1 = map(radians, coord1)
    lon2, lat2 = map(radians, coord2)
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    km = 6367 * c
    return km
