#!/Users/nicolasf/anaconda/bin/python
import os
import sys
import argparse

from collections import OrderedDict as od
import json

import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, '../')

from paleopy import analogs
from paleopy import ensemble
from paleopy import WR
from paleopy.plotting import scalar_plot
from paleopy.plotting import indices

"""
Global function for saving progress
"""


def save_progress(path=None, step=None, value=0):
    progress = od()
    progress['step'] = step
    progress['percentage'] = value
    with open(os.path.join(path, 'output.json'), 'w') as f:
        json.dump(progress, f)


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

parser.add_argument('-v', '--verbose', dest='verbose', type=bool, default=False, \
                    help='Output progress')

"""
goes from argparse Namespace to a dictionnary or key / value arguments
"""

vargs = vars(parser.parse_args())

"""
pop `opath` (the path where the outputs are saved) out of the dictionnary
"""

opath = vargs.pop('opath')

"""
pop `verbose` out of the dictionnary
"""
verbose = vargs.pop('verbose')


"""
instantiates an `ensemble` class, pass the `vargs` dict of keyword arguments to the class
"""

ens = ensemble(**vargs)

# p = proxy(sitename, lon, lat, dpath=dpath, dataset=dataset, variable=variable,
#           season=season, value=value, period=period, climatology=climatology, calc_anoms=calc_anoms, detrend=detrend)

"""
Creates output file array
"""

images = []

if verbose:
    save_progress(opath, 'SST', 0)
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
images.append({'id': 'SST_GL_map_ensemble', 'title' : 'SST', 'filename': 'map_ensemble.png'})

if verbose:
    save_progress(opath, 'UWND at 200hpa', 20)

# ==============================================================================
"""
UWND at 850 and 200 hPa
"""

uwnd = analogs(ens, 'ncep', 'uwnd_200').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath, 'map2_ensemble.png'))
images.append({'id': 'NCEP_GL_map_ensemble_200', 'title' : 'Zonal Wind at 200hPa', 'filename': 'map2_ensemble.png'})

if verbose:
    save_progress(opath, 'UWND at 850hpa', 40)

uwnd = analogs(ens, 'ncep', 'uwnd_850').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath, 'map3_ensemble.png'))
images.append({'id': 'NCEP_GL_map_ensemble_850', 'title' : 'Zonal Wind at 850hPa', 'filename': 'map3_ensemble.png'})

if verbose:
    save_progress(opath, 'Climate Indices', 60)

# ==============================================================================
"""
CLIMATE INDICES
"""

f = indices(ens).plot()

f.savefig(os.path.join(opath, 'indices_ensemble.png'))
images.append({'id': 'indices_ensemble', 'title' : 'Climate Indices', 'filename': 'indices_ensemble.png'})

if verbose:
    save_progress(opath, 'Weather Regimes', 80)
# ==============================================================================
"""
WEATHER REGIMES
"""

w = WR(ens, classification='New Zealand')

w.probs_anomalies(kind='many')

f = w.plot_bar()

f.savefig(os.path.join(opath, 'WR_ensemble.png'))
images.append({'id': 'Ensemble_bar_plot', 'title' : 'Weather Regimes', 'filename': 'WR_ensemble.png'})

if verbose:
    save_progress(opath, 'Complete', 100)


# Save images list to json file
with open(os.path.join(opath, 'images.json'), 'w') as f:
    json.dump(images, f)