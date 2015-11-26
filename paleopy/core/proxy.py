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

import warnings
warnings.filterwarnings(action='ignore', category=FutureWarning)

class proxy:
    """
    base class for a single proxy

    Parameters
    ----------

    sitename : string
            The name of the proxy site, no default
    lon : float
            The longitude (in decimal degrees) of the site
    lat : float
            The latitude (in decimal degrees) of the site
    dataset : string
            The dataset to load to look for analogs, see
            the `datasets.json` file for a list of the available
            datasets (e.g. `ersst`, 'ncep', `gpcp`)
            default is `ersst`
    variable : string
            The variable to extract from the dataset to look for analogs, see
            the `datasets.json` file for a list of variables available for
            each dataset (e.g. for `ncep` dataset: `hgt_1000, hgt_850` etc
            default is `sst`
    season : string
            The season to which the proxy is sensitive, can be:
            - any of `DJF`, `JFM`, `FMA`, ..
            - `Warm Season (Dec. - May)`
            - `Cold Season (Jun. - Nov.)`
            - `Year (Jan. - Dec.)`
            - `Hydro. year (Jul. - Jun.)`
            default is `DJF`
    value : float or integer or string
            The value attached to the proxy
            if string, must be in ['WB', 'B', 'N', 'A', 'WA']
            for Well-Below, Below, etc and will be taken
            as the category of anomalies WRT to present conditions
    period : tuple of integers
            The period from which the analogs can be taken in a
            given dataset. Default is (1979, 2014)
    climatology : tuple of integers
            The climatological period with respect to which the
            anomalies (if `calc_anoms == True`) are calculated
    calc_anoms : boolean
            If `True`, the anomalies are calculated before the analog
            years are determined: The `value` parameter needs to represent
            an *anomaly* with respect to present condition
            If `False`, the analog years are determined from the raw seasonal
            time-series: The `value` parameter needs to represent a raw value
            (i.e. a rainfall amount, a mean temperature)
    detrend : boolean
            If `True` the linear trend is removed from the time-series
            before the analog years are determined, if `False`, nothing is
            done

    Attributes
    ----------
    """

    def __init__(self, sitename=None, lon=None, lat=None, dpath='./jsons', dataset='ersst', variable='sst', season='DJF', value=None, \
                 period=(1979, 2014), climatology=(1981,2010), calc_anoms=True, detrend=True):
        super(proxy, self).__init__()
        if lon < 0:
            lon += 360.
        self.description = 'proxy'
        self.sitename = sitename
        self.coords = (lon, lat)
        self.jsons = dpath
        self.dataset = dataset
        self.variable = variable
        self.season = season
        self.value = value
        self.period = period
        self.climatology = climatology
        self.calc_anoms = calc_anoms
        self.detrend = detrend

    def read_dset_params(self):
        with open(os.path.join(self.jsons, 'datasets.json'), 'r') as f:
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
            self.ts = pd.DataFrame(ts.loc[:,self.variable])
        else:
            lon = dset['longitudes'].data
            lat = dset['latitudes'].data
            lons, lats = np.meshgrid(lon, lat)
            lonf = lons.flatten('F')
            latf = lats.flatten('F')
            ### TODO:
            ### test, then replace do_kdtree with the call to the sel
            ### method of a xray.Dataset with method = 'nearest'
            self.extracted_coords = do_kdtree(lonf, latf, point)
            self.distance_point = haversine(self.extracted_coords, point)
            ts = dset[self.variable].sel(longitudes=self.extracted_coords[0], latitudes=self.extracted_coords[1])
            ts = ts.to_dataframe()
            self.ts = pd.DataFrame(ts.loc[:,self.variable])
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

        # casts that into a pandas DataFrame
        ts_seas.loc[:,'anomalies'] = ts_seas - ts_seas.ix[start_clim:end_clim].mean(0)

        # caculates trend for raw values
        x = ts_seas.index.year
        y = ts_seas.loc[:,self.variable]
        slope, intercept, pval, rval, stderr = linregress(x, y)
        if not self.calc_anoms:
            self.trend_params = {'slope':slope, 'intercept':intercept}
        yhat = slope * x + intercept
        ydetrend = y - yhat
        ts_seas.loc[:,'d_' + self.variable] = (ydetrend + y.mean())

        # calculates the trend for the anomalies
        x = ts_seas.index.year
        y = ts_seas.loc[:,'anomalies']
        slope, intercept, pval, rval, stderr = linregress(x, y)
        if self.calc_anoms:
            self.trend_params = {'slope':slope, 'intercept':intercept}
        yhat = slope * x + intercept
        ydetrend = y - yhat
        ts_seas.loc[:,'d_anomalies'] = ydetrend

        self.ts_seas = ts_seas

    def find_analogs(self):
        val = self.value
        if self.calc_anoms and not self.detrend:
            ts = self.ts_seas.loc[:,'anomalies']
        if self.calc_anoms and self.detrend:
            ts = self.ts_seas.loc[:,'d_anomalies']
        if not self.calc_anoms and self.detrend:
            ts = self.ts_seas.loc[:,'d_' + self.variable]
        if not self.calc_anoms and not self.detrend:
            ts = self.ts_seas.loc[:,self.variable]

        labels=['WB','B','N','A','WA']
        self.ts_seas.loc[:,'cat'], bins = pd.qcut(ts, 5, labels=labels, retbins=True)
        # if value passed is a string, we get from the category
        if isinstance(val, str):
            subset = self.ts_seas[self.ts_seas['cat'] == val]
            self.category = val
        # if not, then we search in the bins
        else:
            bins[0] = -np.inf
            bins[-1] = np.inf
            category = labels[np.searchsorted(bins, val)-1]
            subset = self.ts_seas[self.ts_seas['cat'] == category]
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
            json_path = os.path.join(self.jsons, 'proxies')
            proxy_name = self.sitename.replace(" ","_")
            proxy_name = proxy_name.replace(".","")
            #proxy_name =
            fname = "{}.json".format(self.sitename.replace(" ","_"))
            with open(os.path.join(json_path, fname),'w') as f:
                json.dump(proxy_dict, f)
        self.proxy_dict = proxy_dict

    def plot_season_ts(self, fname=None):

        f, ax = plt.subplots(figsize=(8,5))

        if self.calc_anoms:
            y = self.ts_seas.loc[:,'anomalies']
            vmin = y.min() -0.1 * (np.abs(y.min()))
            vmax = y.max() +0.1 * (np.abs(y.min()))

            yd = self.ts_seas.loc[:,'d_anomalies']

            if self.detrend:
                ya = self.analogs.loc[:,'d_anomalies']
            else:
                ya = self.analogs.loc['anomalies']

        else:

            y = self.ts_seas.loc[:,self.variable]
            vmin = y.min() -0.01 * (np.abs(y.min()))
            vmax = y.max() +0.01 * (np.abs(y.min()))

            yd = self.ts_seas.loc[:,'d_' + self.variable]

            if self.detrend:
                ya = self.analogs.loc[:,'d_' + self.variable]
            else:
                ya = self.analogs.loc[self.variable]

        ax.plot(y.index, y.values, 'steelblue', lw=2, label='ts')
        ax.plot(yd.index, yd.values, color='k', lw=2, label='ts (detrended)')
        ax.plot(ya.index, ya.values, 'ro', label='analog years')
        ax.vlines(ya.index, vmin, vmax, lw=5, alpha=0.4, label="")

        ax.set_xlim(y.index[0] - relativedelta(years=1), y.index[-1] + relativedelta(years=1))

        ax.set_ylim(vmin, vmax)
        ax.set_ylabel(self.variable +": " + self.dset_dict['units'], fontsize=14)

        ax.legend(framealpha=0.4, loc='best')

        [ax.axhline(b, color='magenta', zorder=1, alpha=0.5) for b in self.quintiles[1:-1]]

        ax.grid()

        lyears = ",".join(map(str, self.analogs.index.year.tolist()))

        ax.set_title("Analog seasons for {} {} from {} {}:\n{}".format(self.season, self.sitename, self.dataset,
                                                                    self.variable, lyears, fontsize=14))

        [l.set_fontsize(14) for l in ax.xaxis.get_ticklabels()]
        [l.set_fontsize(14) for l in ax.yaxis.get_ticklabels()]

        if fname:
            f.savefig(fname, dpi=200)
            plt.close(f)

        return f
