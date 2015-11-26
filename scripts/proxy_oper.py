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

parser.add_argument('-d','--dpath', dest='dpath', type=str, default='../jsons', \
help='the path to the data (jsons and associated netcdf and csv files)')

parser.add_argument('-n','--name', dest='sitename', type=str, default='Rarotonga', \
help='the name of the site')

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

parser.add_argument('-an','--calc_anoms', dest='calc_anoms', type=bool, default=True, \
help='True if the anomalies are calculated, False otherwise. Default is True')

parser.add_argument('-dt','--detrend', dest='detrend', type=bool, default=True, \
help='True if the time-series need detrended, False otherwise. Default is True')

args = parser.parse_args()

globals().update(vars(args))

"""
instantiate a proxy class
"""

p = proxy(sitename, lon, lat, dpath=dpath, dataset=dataset, variable=variable,
          season=season, value=value, calc_anoms=calc_anoms, detrend=detrend)


"""
process the proxy
"""


p.extract_ts()
p.calculate_season()
p.find_analogs()
f = p.plot_season_ts()
f.savefig('./essai.png')

"""
instantiate the analog classes with the proxy for each dataset + variable we
want to map
"""

sst = analogs(p, 'ersst', 'sst').composite()

f = scalar_plot(sst, test=0.1, proj='cyl').plot()

f.savefig('./map.png')

uwnd = analogs(p, 'ncep', 'uwnd_200').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig('./map2.png')

uwnd = analogs(p, 'ncep', 'uwnd_850').composite()

f = scalar_plot(uwnd, test=0.05, proj='cyl').plot()

f.savefig('./map3.png')

"""
instantiate an `indices` class and plot
"""

f = indices(p).plot()
f.savefig('./indices.png')
