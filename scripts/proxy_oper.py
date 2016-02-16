#!/Users/nicolasf/anaconda/bin/python
import os
import sys
import argparse

import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, '../')

from paleopy import proxy
from paleopy import analogs
from paleopy import ensemble
from paleopy.plotting import scalar_plot
from paleopy.plotting import indices

"""
parse command line arguments
"""

parser = argparse.ArgumentParser()

parser.add_argument('-dj','--djsons', dest='djsons', type=str, default='../jsons', \
help='the path to the jsons files defining the paths and parameters of each dataset')

parser.add_argument('-pj','--pjsons', dest='pjsons', type=str, default='../jsons/proxies', \
help='the path where to save the individual proxy json files')

parser.add_argument('-jn','--pfname', dest='pfname', type=str, default=None, \
help='the name of the JSON file containing the information for a single proxy')

parser.add_argument('-o','--opath', dest='opath', type=str, default='./outputs', \
help='the path where to save the figures')

parser.add_argument('-n','--name', dest='sitename', type=str, default='Rarotonga', \
help='the name of the site')

parser.add_argument('-t','--type', dest='proxy_type', type=str, default='Coral core', \
help='the type of proxy (coral, Tree-ring, etc)')

parser.add_argument('-lon','--longitude', dest='lon', type=float, default=-159.82, \
help='the longitude (decimal degree) of the proxy site')

parser.add_argument('-lat','--latitude', dest='lat', type=float, default=-21.23, \
help='the latitude (decimal degree) of the proxy site')

parser.add_argument('-dset','--dataset', dest='dataset', type=str, default='ersst', \
help='the dataset to interrogate to draw the analog years')

parser.add_argument('-var','--variable', dest='variable', type=str, default='sst', \
help='the variable in the dataset to interrogate to draw the analog years')

parser.add_argument('-s','--season', dest='season', type=str, default='DJF', \
help='the season to which the proxy is sensitive')

parser.add_argument('-val','--value', dest='value', default=0.6, \
help='the value for the proxy (can be either a float or a string)')

parser.add_argument('-per','--period', dest='period', type=tuple, default=(1979, 2014), \
help='the period from which to draw the analog seasons')

parser.add_argument('-clim','--climatology', dest='climatology', type=tuple, default=(1981, 2010), \
help='the climatological period with respect to which the anomalies are calculated')

parser.add_argument('-an','--calc_anoms', dest='calc_anoms', type=bool, default=True, \
help='True if the anomalies are calculated, False otherwise. Default is True')

parser.add_argument('-dt','--detrend', dest='detrend', type=bool, default=True, \
help='True if the time-series need detrended, False otherwise. Default is True')

# new arguments as from the 27 of January 2016

parser.add_argument('-a','--aspect', dest='aspect', type=float, default=None, \
help='the aspect (in degrees, from 0 to 360)')

parser.add_argument('-e','--elevation', dest='elevation', type=float, default=None, \
help='the elevation (in meters)')

parser.add_argument('-dc','--dating', dest='dating_convention', type=str, default=None, \
help='the dating convention')

parser.add_argument('-cal','--calendar', dest='calendar', type=str, default=None, \
help='the calendar year')

parser.add_argument('-ch','--chronology', dest='chronology', type=str, default=None, \
help='the chronology control (i.e. 14C, Historic, Dendrochronology, etc)')

parser.add_argument('-m','--measurement', dest='measurement', type=str, default=None, \
help='the proxy measurement type (e.g. width for tree rings)')

"""
goes from argparse Namespace to a dictionnary or key / value arguments
"""

vargs = vars(parser.parse_args())

"""
pop `opath` (the path where the outputs are saved) out of the dictionnary
"""

opath = vargs.pop('opath')

"""
instantiates a proxy class, pass the `vargs` dict of keyword arguments to the class
"""

p = proxy(**vargs)

"""
process the proxy
"""

p.extract_ts()
p.calculate_season()
p.find_analogs()
f = p.plot_season_ts()
p.proxy_repr()
f.savefig(os.path.join(opath,'time_series.png'))

"""
instantiate the analog classes with the proxy for each dataset + variable we
want to map
"""

# ==============================================================================
"""
SST
"""

sst = analogs(p, 'ersst', 'sst').composite()

f = scalar_plot(sst, test=0.1, proj='cyl').plot()

f.savefig(os.path.join(opath,'map1_proxy.png'))

# ==============================================================================
"""
UWND at 850 and 200 hPa
"""

uwnd = analogs(p, 'ncep', 'uwnd_200').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath,'map2_proxy.png'))

uwnd = analogs(p, 'ncep', 'uwnd_850').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath,'map3_proxy.png'))

# ==============================================================================
"""
CLIMATE INDICES
"""

f = indices(p).plot()

f.savefig(os.path.join(opath,'indices_proxy.png'))
