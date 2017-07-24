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
from paleopy.plotting import scalar_plot
from paleopy.plotting import indices
from paleopy.utils import save_progress
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
if verbose:
    save_progress(opath, 'Process the proxy', 0)
p.find_analogs()

# 2: save the proxy in the JSON file
p.proxy_repr()

"""
plots the time series for the proxy analogs
"""
if verbose:
    save_progress(opath, 'Time Series', 10)
f = p.plot_season_ts()

f.savefig(os.path.join(opath, 'time_series.png'))
images.append({'id': 'time_series', 'title' : 'Analog Seasons', 'filename': 'time_series.png'})
plt.close(f)

"""
plots the VCSN composite anomalies in every case
"""

if verbose:
    save_progress(opath, 'VCSN', 15)

vcsn = Analogs(p, 'vcsn', 'Rain').composite()
f = scalar_plot(vcsn, test=0.1, proj='cyl', res='h').plot(subplots=False)
f.savefig(os.path.join(opath,'VCSN_rain_proxy.png'))
images.append({'id': 'vcsn_rain', 'title' : 'VCSN seasonal rainfall', 'filename': 'VCSN_rain_proxy.png'})
plt.close(f)

vcsn = Analogs(p, 'vcsn', 'TMean').composite()
f = scalar_plot(vcsn, test=0.1, proj='cyl', res='h', vmin=-1.5, vmax=1.5).plot(subplots=False)
f.savefig(os.path.join(opath,'VCSN_tmean_proxy.png'))
images.append({'id': 'vcsn_tmean', 'title' : 'VCSN seasonal Temperatures', 'filename': 'VCSN_tmean_proxy.png'})
plt.close(f)

vcsn = Analogs(p, 'vcsn', 'Vapour_Pressure').composite()
f = scalar_plot(vcsn, test=0.1, proj='cyl', res='h', vmin=-3, vmax=3).plot(subplots=False)
f.savefig(os.path.join(opath,'VCSN_VP_proxy.png'))
images.append({'id': 'vcsn_VP', 'title' : 'VCSN seasonal Vapour Pressure', 'filename': 'VCSN_VP_proxy.png'})
plt.close(f)


# ==============================================================================
"""
Sea Surface Temperatures, global
"""
if verbose:
    save_progress(opath, 'SST', 20)
sst = Analogs(p, 'ersst', 'sst').composite()
f = scalar_plot(sst, test=0.1, proj='cyl', vmin=-1.5, vmax=1.5).plot()
f.savefig(os.path.join(opath,'SST_proxy.png'))
images.append({'id': 'sst', 'title' : 'Sea Surface Temperature', 'filename': 'SST_proxy.png'})
plt.close(f)


"""
HGT at 850 hPa, global
"""
if verbose:
    save_progress(opath, 'HGT at 850 hPa', 30)
hgt = Analogs(p, 'ncep', 'hgt_850').composite()
f = scalar_plot(hgt, test=0.05, proj='cyl', vmin=-50.0, vmax=50.0).plot()
f.savefig(os.path.join(opath,'hgt_850_proxy.png'))
images.append({'id': 'hgt_850', 'title' : 'Geopotential at 850 hPa', 'filename': 'hgt_850_proxy.png'})
plt.close(f)

"""
HGT at 1000 hPa, global
"""
if verbose:
    save_progress(opath, 'HGT at 1000 hPa', 32)
hgt = Analogs(p, 'ncep', 'hgt_1000').composite()
f = scalar_plot(hgt, test=0.05, proj='cyl', vmin=-25.0, vmax=25.0).plot()
f.savefig(os.path.join(opath,'hgt_1000_proxy.png'))
images.append({'id': 'hgt_1000', 'title' : 'Geopotential at 1000 hPa', 'filename': 'hgt_1000_proxy.png'})
plt.close(f)

"""
HGT at 1000 hPa, polar
"""
if verbose:
    save_progress(opath, 'HGT at 1000 hPa', 34)
hgt = Analogs(p, 'ncep', 'hgt_1000').composite()
f = scalar_plot(hgt, test=0.05, proj='spstere', vmin=-25.0, vmax=25.0).plot()
f.savefig(os.path.join(opath,'hgt_1000_proxy_polar.png'))
images.append({'id': 'hgt_1000', 'title' : 'Geopotential at 1000 hPa', 'filename': 'hgt_1000_proxy_polar.png'})
plt.close(f)

"""
HGT at 1000 hPa, NZ domain, composite
"""

if verbose:
    save_progress(opath, 'HGT 1000 global', 36)
f = scalar_plot(hgt, test=0.1, proj='cyl', domain=[165, 180, -50., -30], res='h', vmin=-25.0, vmax=25.0).plot(subplots=False)
f.savefig(os.path.join(opath,'hgt_1000_proxy_NZ.png'))
images.append({'id': 'hgt_1000_NZ', 'title' : 'Geopotential at 1000 hPa, NZ domain', 'filename': 'hgt_1000_proxy_NZ.png'})
plt.close(f)

"""
HGT at 1000 hPa, NZ domain, one map per year
"""
if verbose:
    save_progress(opath, 'HGT 1000 NZ domain Analogs', 38)
f = scalar_plot(hgt, test=0.1, proj='cyl', domain=[165, 180, -50., -30], res='h', vmin=-25.0, vmax=25.0).plot(subplots=True)
f.savefig(os.path.join(opath,'hgt_1000_proxy_NZ_years.png'))
images.append({'id': 'hgt_1000_NZ_samples', 'title' : 'Geopotential at 1000 hPa, NZ domain, Analog years', 'filename': 'hgt_1000_proxy_NZ_years.png'})
plt.close(f)

"""
UWND at 200 hPa
"""
if verbose:
    save_progress(opath, 'UWND at 200hpa', 40)
uwnd = Analogs(p, 'ncep', 'uwnd_200').composite()
f = scalar_plot(uwnd, test=0.05, proj='cyl', vmin=-10.0, vmax=10.0).plot()
f.savefig(os.path.join(opath, 'map2_proxy.png'))
images.append({'id': 'uwnd_200', 'title' : 'Zonal Wind at 200 hPa', 'filename': 'map2_proxy.png'})
plt.close(f)

"""
UWND at 850 hPa
"""
if verbose:
    save_progress(opath, 'UWND at 850hpa', 43)
uwnd = Analogs(p, 'ncep', 'uwnd_850').composite()
f = scalar_plot(uwnd, test=0.05, proj='cyl', vmin=-3.0, vmax=3.0).plot()
f.savefig(os.path.join(opath, 'map3_proxy.png'))
images.append({'id': 'uwnd_850', 'title' : 'Zonal Wind at 850 hPa', 'filename': 'map3_proxy.png'})
plt.close(f)

"""
UWND at 1000 hPa
"""
if verbose:
    save_progress(opath, 'UWND at 1000hpa', 46)
uwnd = Analogs(p, 'ncep', 'uwnd_1000').composite()
f = scalar_plot(uwnd, test=0.05, proj='cyl', vmin=-3.0, vmax=3.0).plot()
f.savefig(os.path.join(opath, 'map4_proxy.png'))
images.append({'id': 'uwnd_1000', 'title' : 'Zonal Wind at 1000 hPa', 'filename': 'map4_proxy.png'})
plt.close(f)

"""
HGT and WIND at 850 hPa
"""
if verbose:
    save_progress(opath, 'hgt and wind at 850hpa', 46)

hgt = Analogs(p, 'ncep', 'hgt_850').composite()
uwnd = Analogs(p, 'ncep', 'uwnd_850').composite()
vwnd = Analogs(p, 'ncep', 'vwnd_850').composite()

f = vector_plot(uwnd, vwnd, hgt).plot(stepp=3, res='c')

f.savefig(os.path.join(opath, 'hgt_wind_850.png'))
images.append({'id': 'hgt_wind_850', 'title' : 'geopotential and wind at 850 hPa', 'filename': 'hgt_wind_850_proxy.png'})
plt.close(f)

"""
HGT and WIND at 200 hPa
"""
if verbose:
    save_progress(opath, 'hgt and wind at 850hpa', 46)

hgt = Analogs(p, 'ncep', 'hgt_200').composite()
uwnd = Analogs(p, 'ncep', 'uwnd_200').composite()
vwnd = Analogs(p, 'ncep', 'vwnd_200').composite()

f = vector_plot(uwnd, vwnd, hgt).plot(stepp=3, res='c')

f.savefig(os.path.join(opath, 'hgt_wind_200.png'))
images.append({'id': 'hgt_wind_200', 'title' : 'geopotential and wind at 200 hPa', 'filename': 'hgt_wind_200_proxy.png'})
plt.close(f)

# ==============================================================================
"""
CLIMATE INDICES
"""
if verbose:
    save_progress(opath, 'climate indices', 50)
f = indices(p).plot()
f.savefig(os.path.join(opath, 'indices_proxy.png'))
images.append({'id': 'indices_proxy', 'title' : 'Climate Indices', 'filename': 'indices_proxy.png'})
plt.close(f)

"""
NZ Weather Regimes
"""
if verbose:
    save_progress(opath, 'NZ weather regimes', 60)
w = WR(p, classification='New Zealand')
f = w.plot_bar(sig=1)
f.savefig(os.path.join(opath, 'NZ_regimes_proxy.png'))
images.append({'id': 'NZ_regimes_proxy', 'title' : 'NZ weather regimes (Kidson Types)', 'filename': 'NZ_regimes_proxy.png'})
plt.close(f)

"""
SW Pacific Weather Regimes
"""
if verbose:
    save_progress(opath, 'SW Pacific weather regimes', 70)
w = WR(p, classification='SW Pacific')
f = w.plot_bar(sig=1)
f.savefig(os.path.join(opath, 'SWPac_regimes_proxy.png'))
images.append({'id': 'SWPac_regimes_proxy', 'title' : 'Southwest Pacific weather regimes', 'filename': 'SWPac_regimes_proxy.png'})
plt.close(f)


if verbose:
    save_progress(opath, 'Complete', 100)

# Save images list to json file
with open(os.path.join(opath, 'images.json'), 'w') as f:
    json.dump(images, f)
