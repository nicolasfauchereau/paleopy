#!/Users/nicolasf/anaconda/bin/python
import os
import sys
import argparse

import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, '../')

from paleopy import analogs
from paleopy import ensemble
from paleopy import WR
from paleopy.plotting import scalar_plot
from paleopy.plotting import indices

"""
parse command line arguments
"""

parser = argparse.ArgumentParser()

parser.add_argument('-dj','--djsons', dest='djsons', type=str, default='../jsons', \
help='the path to the jsons files defining the paths and parameters of each dataset')

parser.add_argument('-j','--pjsons', dest='pjsons', type=str, default='../jsons/proxies', \
help='the directory containing the proxy json files')

parser.add_argument('-o','--opath', dest='opath', type=str, default='./outputs', \
help='the directory in which to save the figures')

parser.add_argument('-s','--season', dest='season', type=str, default='DJF', \
help='the season to consider: will be checked against the individual proxies seasons')


"""
goes from argparse Namespace to a dictionnary or key / value arguments
"""

vargs = vars(parser.parse_args())

opath = vargs.pop('opath')

"""
instantiates an `ensemble` class, pass the `vargs` dict of keyword arguments to the class
"""

ens = ensemble(**vargs)

# p = proxy(sitename, lon, lat, dpath=dpath, dataset=dataset, variable=variable,
#           season=season, value=value, period=period, climatology=climatology, calc_anoms=calc_anoms, detrend=detrend)


"""
instantiate the analog classes with the proxy for each dataset + variable we
want to map
"""

# ==============================================================================
"""
SST
"""

sst = analogs(ens, 'ersst', 'sst').composite()

f = scalar_plot(sst, test=0.1, proj='cyl').plot()

f.savefig(os.path.join(opath, 'map_ensemble.png'))

# ==============================================================================
"""
UWND at 850 and 200 hPa
"""

uwnd = analogs(ens, 'ncep', 'uwnd_200').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath, 'map2_ensemble.png'))

uwnd = analogs(ens, 'ncep', 'uwnd_850').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath, 'map3_ensemble.png'))

# ==============================================================================
"""
CLIMATE INDICES
"""

f = indices(ens).plot()

f.savefig(os.path.join(opath, 'indices_ensemble.png'))

# ==============================================================================
"""
WEATHER REGIMES
"""

w = WR(ens, classification='Kidson Types')

w.probs_anomalies(kind='many')

f = w.plot_bar()

f.savefig(os.path.join(opath, 'WR_ensemble.png'))
