# Python packages imports
import os
from collections import OrderedDict as od
import numpy as np
from numpy import ma
import pandas as pd
import matplotlib.pyplot as plt
from dateutil.relativedelta import relativedelta
import json
import xray
from scipy.stats import linregress

# relative imports
from ..utils.do_kdtree import do_kdtree
from ..utils.haversine import haversine
from ..utils.pprint_od import pprint_od


class proxy:
    'base class for a single proxy'

    def __init__(self, sitename, lon, lat, dataset, variable, season, value=None, \
                 period=(1979, 2014), climatology=(1981,2010), calc_anoms=True, detrend=True):
        super(proxy, self).__init__()
        if lon < 0:
            lon += 360.
        self.description = 'proxy'
        self.sitename = sitename
        self.coords = (lon, lat)
        self.dataset = dataset
        self.variable = variable
        self.season = season
        self.value = value
        self.period = period
        self.climatology = climatology
        self.calc_anoms = calc_anoms
        self.detrend = detrend

    def read_dset_params(self, dpath='./jsons'):
        with open(os.path.join(dpath, 'datasets.json'), 'r') as f:
            dset_dict = json.loads(f.read())
        self.dset_dict = dset_dict[self.dataset][self.variable]

    def check_domain(self):
        if not(hasattr(self, 'dset_dict')):
            self.read_dset_params()
        domain = self.dset_dict['domain']
        lond = self.coords[0]
        latd = self.coords[1]
        # uncomment if one wants to print
        # print("\nLON: {:4.2f}, LAT: {:4.2f}".format(lond, latd))
        if ( (lond <= domain[0]) \
            | (lond >= domain[1]) \
            | (latd <= domain[2]) \
            | (latd >= domain[3]) ):
            print("""
            ERROR! coordinates of the proxy fall outside
            of the bounds of the domain for dataset {}
            """.format(self.dataset))
            raise Exception("DOMAIN ERROR")

    def extract_ts(self):
        # checks the domain first
        self.check_domain()
        # if all good, we proceed
        fname = self.dset_dict['path']
        point = self.coords
        start = str(self.period[0])
        end = str(self.period[1])
        dset = xray.open_dataset(fname)
        dset = dset.sel(time=slice(start, end))
        # test is a mask is present, if so assume we have a grid
        # then meshgrid, mask, flatten, etc
        if 'mask' in dset.data_vars:
            mask = dset['mask'].data
            mask = mask.astype(np.bool)
            lon = dset['longitudes'].data
            lat = dset['latitudes'].data
            lons, lats = np.meshgrid(lon, lat)
            lons = ma.masked_array(lons, mask)
            lats = ma.masked_array(lats, mask)
            lonf = lons.flatten('F').compressed()
            latf = lats.flatten('F').compressed()
            self.extracted_coords = do_kdtree(lonf, latf, point)
            self.distance_point = haversine(self.extracted_coords, point)
            ts = dset[self.variable].sel(longitudes=self.extracted_coords[0], latitudes=self.extracted_coords[1])
            ts = ts.to_dataframe()
            self.ts = ts.loc[:,self.variable]
        else:
            lon = dset['longitudes'].data
            lat = dset['latitudes'].data
            lons, lats = np.meshgrid(lon, lat)
            lonf = lons.flatten('F')
            latf = lats.flatten('F')
            self.extracted_coords = do_kdtree(lonf, latf, point)
            self.distance_point = haversine(self.extracted_coords, point)
            ts = dset[self.variable].sel(longitudes=self.extracted_coords[0], latitudes=self.extracted_coords[1])
            ts = ts.to_dataframe()
            self.ts = ts.loc[:,self.variable]
        dset.close()

    def calculate_season(self):
        season = self.season
        start_clim = str(self.climatology[0])
        end_clim = str(self.climatology[1])
        # seasons parameters is a dictionnary with:
        # key = the season string ('DJF', 'JJA')
        # value =  a tuple (length of the season, month of the end of the season)
        seasons_params = {}
        seasons_params['DJF'] = (3,2)
        seasons_params['JFM'] = (3,3)
        seasons_params['FMA'] = (3,4)
        seasons_params['MAM'] = (3,5)
        seasons_params['AMJ'] = (3,6)
        seasons_params['MJJ'] = (3,7)
        seasons_params['JJA'] = (3,8)
        seasons_params['JAS'] = (3,9)
        seasons_params['ASO'] = (3,10)
        seasons_params['SON'] = (3,11)
        seasons_params['OND'] = (3,12)
        seasons_params['NDJ'] = (3,1)
        seasons_params['Warm Season (Dec. - May)'] = (6, 5)
        seasons_params['Cold Season (Jun. - Nov.)'] = (6, 11)
        #seasons_params['Annual'] = (12, 12)
        seasons_params['Year (Jan. - Dec.)'] = (12, 12)
        seasons_params['Hydro. year (Jul. - Jun.)'] = (12, 6)

        # if the variable is rainfall, we calculate rolling sum
        if self.dset_dict['units'] in ['mm']:
            ts_seas = pd.rolling_sum(self.ts, seasons_params[season][0])
        # else we calculate the rolling mean (average)
        else:
            ts_seas = pd.rolling_mean(self.ts, seasons_params[season][0])

        ts_seas = ts_seas[ts_seas.index.month == seasons_params[season][1]]

        # drop the missing values coming from the rolling average
        ts_seas.dropna(inplace=True)

        self.ts_seas = ts_seas

        if self.calc_anoms:
            ts_seas = ts_seas - ts_seas.ix[start_clim:end_clim].mean(0)
            self.ts_seas = ts_seas
            # we calculate the trend but don't subtract it
            x = self.ts_seas.index.year
            y = self.ts_seas.values
            slope, intercept, pval, rval, stderr = linregress(x, y)
            yhat = slope * x + intercept
            self.trend_params = {'slope':slope, 'intercept':intercept, 'pval':pval, 'rval':rval, 'sterr':stderr}
            self.trend_line = pd.DataFrame(yhat, index=self.ts_seas.index, columns=[self.variable+'_trend'])

        if self.detrend:
            # we calculate the trend AND subtract it
            x = self.ts_seas.index.year
            y = self.ts_seas.values
            slope, intercept, pval, rval, stderr = linregress(x, y)
            yhat = slope * x + intercept
            ydetrend = y - yhat
            ts_seas = pd.DataFrame(ydetrend, index=self.ts_seas.index, columns=[self.variable])
            self.ts_seas = ts_seas
            self.trend_params = {'slope':slope, 'intercept':intercept, 'pval':pval, 'rval':rval, 'sterr':stderr}
            self.trend_line = pd.DataFrame(yhat, index=self.ts_seas.index, columns=[self.variable+'_trend'])

        else:
            # we calculate the trend on the raw values, and don't subtract it
            x = self.ts_seas.index.year
            y = self.ts_seas.values
            slope, intercept, pval, rval, stderr = linregress(x, y)
            yhat = slope * x + intercept
            self.trend_params = {'slope':slope, 'intercept':intercept, 'pval':pval, 'rval':rval, 'sterr':stderr}
            self.trend_line = pd.DataFrame(yhat, index=self.ts_seas.index, columns=[self.variable+'_trend'])

        self.ts_seas = pd.DataFrame(self.ts_seas.values, index=self.ts_seas.index, columns=[self.variable])

    def find_analogs(self):
        val = self.value
        ts_seas = self.ts_seas
        labels=['WB','B','N','A','WA']
        ts_seas.loc[:,'cat'], bins = pd.qcut(ts_seas[self.variable], 5, labels=labels, retbins=True)
        # if value passed is a string, we get from the category
        if isinstance(val, str):
            subset = ts_seas[ts_seas['cat'] == val]
            self.category = val
        # if not, then we search in the bins
        else:
            bins[0] = -np.inf
            bins[-1] = np.inf
            category = labels[np.searchsorted(bins, val)-1]
            subset = ts_seas[ts_seas['cat'] == category]
            self.category = category
        self.analogs = subset
        self.analog_years = subset.index.year
        self.quintiles = bins

    def proxy_repr(self, pprint=False, outfile=True, json_path='./jsons/proxies'):
        """
        proxy_dict is an OrderedDict
        """
        proxy_dict = od()
        proxy_dict['sitename'] = self.sitename
        proxy_dict['coords'] = self.coords
        proxy_dict['season'] = self.season
        proxy_dict['dataset'] = self.dataset
        proxy_dict['variable'] = self.variable
        proxy_dict['calc_anoms'] = self.calc_anoms
        proxy_dict['detrend'] = self.detrend
        proxy_dict['value'] = self.value
        proxy_dict['climatology'] = self.climatology
        proxy_dict['period'] = self.period
        proxy_dict['extracted_coords'] = self.extracted_coords.tolist()
        proxy_dict['distance_point'] = self.distance_point
        proxy_dict['trend_params'] = self.trend_params
        proxy_dict['category'] = self.category
        proxy_dict['analog_years'] = self.analog_years.tolist()

        if pprint:
            pprint_od(proxy_dict)

        if outfile:
            proxy_name = self.sitename.replace(" ","_")
            proxy_name = proxy_name.replace(".","")
            #proxy_name =
            fname = "{}.json".format(self.sitename.replace(" ","_"))
            with open(os.path.join(json_path, fname),'w') as f:
                json.dump(proxy_dict, f)
        self.proxy_dict = proxy_dict

    def plot_season_ts(self, fname=None):

        y = self.ts_seas
        analogs = self.analogs
        analog_years = self.analog_years
        variable = self.variable
        vmin = y[variable].min() + (0.1 * y[variable].min())
        vmax = y[variable].max() + (0.1 * y[variable].max())
        trend = self.trend_line
        btrend = self.detrend
        units = self.dset_dict['units']

        f, ax = plt.subplots(figsize=(8,5))

        # if the user chose to detrend, we show the non-detrended time-series
        ax.plot(trend.index, trend.values, '0.4', label='trend')

        if btrend:
            ax.plot(y.index, y[variable] + trend.values.flatten(), 'steelblue', lw=2, label='ts')

        ax.plot(y.index, y[variable], color='k', lw=2, label='ts (detrended)')

        ax.plot(analogs.index, analogs[variable], 'ro')

        ax.set_ylim(vmin, vmax)

        ax.set_ylabel(variable +": " + units, fontsize=14)

        ax.set_xlim(trend.index[0] - relativedelta(years=1), trend.index[-1] + relativedelta(years=1))

        ax.set_xticks(trend.index[np.arange(1,len(y), 4)])

        ax.legend(framealpha=0.4, loc='best')

        ax.vlines(analogs.index, vmin, vmax, lw=5, alpha=0.4)

        [ax.axhline(b, color='magenta', zorder=1, alpha=0.5) for b in self.quintiles[1:-1]]

        ax.grid()

        lyears = ",".join(map(str, analogs.index.year.tolist()))

        ax.set_title("Analog seasons for {} {} from {} {}:\n{}".format('DJF', self.sitename, self.dataset,
                                                                    variable, lyears, fontsize=14))

        [l.set_fontsize(14) for l in ax.xaxis.get_ticklabels()]
        [l.set_fontsize(14) for l in ax.yaxis.get_ticklabels()]

        if fname:
            f.savefig(fname, dpi=200)
            plt.close(f)
