#!/Users/nicolasf/anaconda/anaconda/bin/python
import os
import sys
import argparse
import json

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

sys.path.insert(0, '../')

from paleopy import proxy
from paleopy import analogs
from paleopy import ensemble
from paleopy.plotting import scalar_plot
from paleopy.plotting import indices
from paleopy import WR


"""
import the little progress file indicator, really
useful only in PICT
"""

from paleopy.utils import save_progress

"""
parse command line arguments
"""

parser = argparse.ArgumentParser()

parser.add_argument('-dj','--djsons', dest='djsons', type=str, default='../jsons', \
help='the path to the jsons files defining the paths and parameters of each dataset')

parser.add_argument('-pj','--pjsons', dest='pjsons', type=str, default='../jsons/proxies', \
help='the path where to save the individual proxy json files')

parser.add_argument('-pn','--pfname', dest='pfname', type=str, default=None, \
help='the name of the JSON file containing the information for a single proxy')

parser.add_argument('-o','--opath', dest='opath', type=str, default='./outputs', \
help='the path where to save the figures, tables and csv files')

parser.add_argument('-n','--name', dest='sitename', type=str, default='Rarotonga', \
help='the name of the site')

parser.add_argument('-t','--ptype', dest='proxy_type', type=str, default='Coral core', \
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
help="""the value for the proxy: can be either a float or a string, if a string, must be in
['WB','B','N','A','WA'] and the `qualitative` flag must be set to True""")

parser.add_argument('-q','--qualitative', dest='qualitative', type=int, default=0, \
help='a flag indicating whether the value passed (see above) is qualitative or not, default to 0 (False): \
i.e. interpret the value as a float')

parser.add_argument('-per','--period', dest='period', type=str, default="1979-2014", \
help='the period from which to draw the analog seasons')

parser.add_argument('-clim','--climatology', dest='climatology', type=str, default="1981-2010", \
help='the climatological period with respect to which the anomalies are calculated')

parser.add_argument('-an','--calc_anoms', dest='calc_anoms', type=int, default=1, \
help='True if the anomalies are calculated, False otherwise. Default is 1 (True)')

parser.add_argument('-dt','--detrend', dest='detrend', type=int, default=1, \
help='True if the time-series need detrended, False otherwise. Default is 1 (True)')

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

parser.add_argument('-v', '--verbose', dest='verbose', type=int, default=0,
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

verbose = bool(vargs.pop('verbose'))

"""
instantiates a proxy class, pass the `vargs` dict of keyword arguments to the class
"""

p = proxy(**vargs)

"""
initialise output file list
"""

images = []


"""
process the proxy
"""

# 1: find the analog seasons
p.find_analogs()

# 2: plot the time-series of seasonal values with the analog years and save
f = p.plot()

f.savefig(os.path.join(opath, 'time_series.png'))

images.append({'id': 'time_series', 'title' : 'Analog Seasons', 'filename': 'time_series.png'})

plt.close(f)

# 3: save the proxy in the JSON file
p.proxy_repr()

"""
instantiate the analog classes with the proxy for each dataset + variable we
want to map
"""

"""
if the attached dataset is the VCSN dataset, we plot the corresponding composite
anomalies for the variable the proxy is sensitive to
"""

dataset = p.dataset
variable = p.variable

compos = analogs(p, dataset, variable).composite()

if dataset in ['vcsn']:
    f = scalar_plot(compos, test=0.1, proj='cyl', res='h').plot(subplots=False)
else:
    f = scalar_plot(compos, test=0.1, proj='cyl', res='i').plot(subplots=False)

images.append({'id': '{}_{}'.format(dataset, variable), 'title' : '{} seasonal {}'.format(dataset, variable), 'filename': '{}_{}_proxy.png'.format(dataset, variable)})

plt.close(f)

data_var = {}
data_var['ncep'] = ['hgt_1000','hgt_850']

# ==============================================================================
"""
Sea Surface Temperatures, global
"""

sst = analogs(p, 'ersst', 'sst').composite()

f = scalar_plot(sst, test=0.1, proj='cyl').plot()

f.savefig(os.path.join(opath,'SST_proxy.png'))

images.append({'id': 'sst', 'title' : 'Sea Surface Temperature', 'filename': 'SST_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'SST global', 10)

"""
HGT at 850 hPa, global
"""

hgt = analogs(p, 'ncep', 'hgt_850').composite()

f = scalar_plot(hgt, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath,'hgt_850_proxy.png'))

images.append({'id': 'hgt_850', 'title' : 'Geopotential at 850 hPa', 'filename': 'HGT_850_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'HGT 850 global', 20)

"""
HGT at 1000 hPa, global
"""

hgt = analogs(p, 'ncep', 'hgt_1000').composite()

f = scalar_plot(hgt, test=0.05, proj='cyl').plot()

f.savefig(os.path.join(opath,'hgt_1000_proxy.png'))

images.append({'id': 'hgt_1000', 'title' : 'Geopotential at 1000 hPa', 'filename': 'HGT_1000_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'HGT 1000 global', 30)

"""
HGT at 1000 hPa, NZ domain, composite
"""

f = scalar_plot(hgt, test=0.1, proj='cyl', domain=[165, 180, -50., -30], res='h').plot(subplots=False)

f.savefig(os.path.join(opath,'hgt_1000_proxy_NZ.png'))

images.append({'id': 'hgt_1000_NZ', 'title' : 'Geopotential at 1000 hPa, NZ domain', 'filename': 'HGT_1000_NZ_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'HGT 1000 NZ domain composite', 40)

"""
HGT at 1000 hPa, NZ domain, one map per year
"""

f = scalar_plot(hgt, test=0.1, proj='cyl', domain=[165, 180, -50., -30], res='h').plot(subplots=True)

f.savefig(os.path.join(opath,'hgt_1000_proxy_NZ_years.png'))

images.append({'id': 'hgt_1000_NZ_samples', 'title' : 'Geopotential at 1000 hPa, NZ domain, analog years', 'filename': 'HGT_1000_NZ_sample_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'HGT 1000 NZ domain analogs', 50)

# ==============================================================================
"""
CLIMATE INDICES
"""

f = indices(p).plot()

f.savefig(os.path.join(opath, 'indices_proxy.png'))

images.append({'id': 'indices_proxy', 'title' : 'Climate Indices', 'filename': 'indices_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'climate indices', 60)


"""
Weather Regimes
"""

w = WR(p, classification='New Zealand')

f = w.plot_bar(sig=1)

f.savefig(os.path.join(opath, 'NZ_regimes_proxy.png'))

images.append({'id': 'NZ_regimes_proxy', 'title' : 'NZ weather regimes (Kidson Types)', 'filename': 'NZ_regimes_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'NZ weather regimes', 70)

w = WR(p, classification='SW Pacific')

f = w.plot_bar(sig=1)

f.savefig(os.path.join(opath, 'SWPac_regimes_proxy.png'))

images.append({'id': 'SWPac_regimes_proxy', 'title' : 'Southwest Pacific weather regimes', 'filename': 'SWPac_regimes_proxy.png'})

plt.close(f)

if verbose:
    save_progress(opath, 'SW Pacific weather regimes', 80)

if verbose:
    save_progress(opath, 'Complete', 100)

# Save images list to json file
with open(os.path.join(opath, 'images.json'), 'w') as f:
    json.dump(images, f)
