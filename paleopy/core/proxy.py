# Python packages imports
import os
from collections import OrderedDict as od
import numpy as np
from numpy import ma
import pandas as pd
import matplotlib.pyplot as plt
import json

try:
    import xarray as xray
except:
    try:
        import xray
    except ImportError:
        print('cannot import xarray or xray')

from scipy.stats import linregress

# relative imports
from ..utils import do_kdtree
from ..utils import haversine
from ..utils import pprint_od
from ..utils import seasons_params
from ..plotting import plot_season_ts

import warnings
warnings.filterwarnings(action='ignore', category=FutureWarning)

class proxy:
    """
    base class for a single proxy

    Parameters
    ----------

    sitename : string
            The name of the proxy site
            user-defined, no default

    proxy_type : string
            the type of proxy
            can be e.g.:
            "Tree-ring"
            "Speleotheme"
            "Coral core"
            user-defined, no default

    lon : float
            The longitude (in decimal degrees) of the site
            user-defined, no default

    lat : float
            The latitude (in decimal degrees) of the site
            user-defined, no default

    djsons : string
            The path to the json files defining the paths
            and parameters arrached to each dataset + variable
            defined by the frontend in PICT, default is ./jsons

    pjsons : string
            The path where to save the individual proxy json files
            defined by the frontend in PICT, default is ./jsons/proxies

    pfname : string
            the name of the JSON file containing the information
            for a single proxy, no default

    dataset : string
            The dataset to load to look for analogs, see
            the `datasets.json` file for a list of the available
            datasets (e.g. `ersst`, 'ncep', `gpcp`)
            user-defined, default is `ersst`

    variable : string
            The variable to extract from the dataset to look for analogs, see
            the `datasets.json` file for a list of variables available for
            each dataset (e.g. for `ncep` dataset: `hgt_1000, hgt_850` etc
            user-defined, default is `sst`

    season : string
            The season to which the proxy is sensitive, can be:
            - any of `DJF`, `JFM`, `FMA`, ..
            - `Warm Season (Dec. - May)`
            - `Cold Season (Jun. - Nov.)`
            - `Year (Jan. - Dec.)`
            - `Hydro. year (Jul. - Jun.)`
            user-defined, default is `DJF`

    value : float or integer or string
            The value attached to the proxy
            if string, must be in ['WB', 'B', 'N', 'A', 'WA']
            for Well-Below, Below, etc and will be taken
            as the category of anomalies WRT to present conditions
            user-defined, no default

    period : string
            The period from which the analogs can be taken in a
            given dataset. Default is "1979-2014"
            user-defined, default is full period for the dataset
            interogated

    climatology : string
            The climatological period with respect to which the
            anomalies (if `calc_anoms == True`) are calculated
            user-defined, default is "1981-2010"

    calc_anoms : boolean
            If `True`, the anomalies are calculated before the analog
            years are determined: The `value` parameter needs to represent
            an *anomaly* with respect to present condition
            If `False`, the analog years are determined from the raw seasonal
            time-series: The `value` parameter needs to represent a raw value
            (i.e. a rainfall amount, a mean temperature)
            user-defined, default is True

    detrend : boolean
            If `True` the linear trend is removed from the time-series
            before the analog years are determined, if `False`, nothing is
            done
            user-defined, default is True

    method : string
            can be either 'closest 8' or 'quintiles'
            to specicify the method employed to choose
            the analog seasons

    aspect : float
            The aspect of the proxy site (in degrees from 0 to 360)
            user-defined, no default

    elevation : float
            The elevation of the proxy site (in meters)
            user-defined, no default

    dating_convention : string
            the dating convention
            user-defined, no default

    calendar : string
            calendar year
            user-defined, no default

    chronology : string
            the chronology control (i.e. 14C, Historic, Dendrochronology, etc)
            user-defined, no default

    measurement : string
            the proxy measurement type (e.g. width for tree rings)
            user-defined, no default

    Attributes
    ----------
    """

    def __init__(self, sitename=None, proxy_type=None, lon=None, lat=None, \
    aspect=None, elevation=None, dating_convention=None, calendar=None, \
    chronology=None, measurement=None, djsons='./jsons', pjsons='./jsons/proxies', \
    pfname=None, dataset='ersst', variable='sst', season='DJF', value=None, \
    qualitative=0, period="1979-2014", climatology="1981-2010", \
    calc_anoms=1, detrend=1, method='quintiles'):

        super(proxy, self).__init__()
        if lon < 0:
            lon += 360.
        self.description = 'proxy'
        self.sitename = sitename
        self.proxy_type = proxy_type
        self.measurement = measurement
        self.coords = (lon, lat)
        self.aspect = aspect
        self.elevation = elevation
        self.djsons = djsons
        self.pjsons = pjsons
        self.pfname = pfname
        self.dataset = dataset
        self.variable = variable
        self.dating_convention = dating_convention
        self.chronology = chronology
        self.calendar = calendar
        self.season = season
        self.value = value
        self.qualitative = bool(qualitative)
        self.period = tuple(map(int,period.split("-"))) # to correct the type
        self.climatology = tuple(map(int,climatology.split("-"))) # to correct the type
        self.calc_anoms = bool(calc_anoms)
        self.detrend = bool(detrend)
        self.method = method

    def read_dset_params(self):
        """
        reads the `datasets.json` file and loads the dictionnary
        containing all the parameters for this dataset
        """
        with open(os.path.join(self.djsons, 'datasets.json'), 'r') as f:
            dset_dict = json.loads(f.read())
        self.dset_dict = dset_dict[self.dataset][self.variable]

    def check_domain(self):
        """
        checks if the domain that is passed
        is compatible with the domain of the
        dataset
        """
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

    def _calc_weights(self, df):
        """
        calculate the weights for compositing
        """
        tmp_df = df.copy(deep=True)
        # print(self.value)
        # print(type(self.value))
        # print(type(tmp_df.iloc[:,0].values))
        # print(np.abs(self.value - tmp_df.iloc[:,0].values.flatten()))
        weights = abs(self.value - tmp_df.iloc[:,0]) / sum(abs(self.value - tmp_df.iloc[:,0]))
        tmp_df.loc[:,'weights'] = (1 - weights) / (1 - weights).sum()
        return tmp_df

    def extract_ts(self):
        """
        extract the time-series for the closest grid-point to
        the passed proxy coordinates
        """
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

        return self

    def calculate_season(self):
        """
        calculates the seasonal values, can be either raw
        or anomalies depending on the parameter (boolean) `calc_anoms`
        passed when instantiating the `proxy` class
        """
        season = self.season
        start_clim = str(self.climatology[0])
        end_clim = str(self.climatology[1])
        # seasons parameters is a dictionnary with:
        # key = the season string ('DJF', 'JJA')
        # value =  a tuple (length of the season, month of the end of the season)
        self.seasons_params = seasons_params()

        if not(hasattr(self, 'ts')):
            self.extract_ts()

        # if the variable is rainfall, we calculate rolling sum
        if self.dset_dict['units'] in ['mm']:
            # test which version of pandas we are using
            if pd.__version__ >= '0.18':
                ts_seas = self.ts.rolling(window=self.seasons_params[season][0]).sum()
            else:
                ts_seas = pd.rolling_sum(self.ts, self.seasons_params[season][0])
        # else we calculate the rolling mean (average)
        else:
            if pd.__version__ >= '0.18':
                ts_seas = self.ts.rolling(window=self.seasons_params[season][0]).mean()
            else:
                ts_seas = pd.rolling_mean(self.ts, self.seasons_params[season][0])

        ts_seas = ts_seas[ts_seas.index.month == self.seasons_params[season][1]]

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

        # TODO: extract the period of interest only, which should
        # be help in self.period
        self.ts_seas = ts_seas.loc[str(self.period[0]):str(self.period[1]),:]

        return self

    def find_analogs(self):
        """
        find the analog seasons

        return:

        self.analogs : a pandas DataFrame
        self.analog_years : a list with the analog years
        self.quintiles : the bins for the quintiles used
        """

        if not(hasattr(self, 'ts_seas')):
             self.calculate_season()

        if self.calc_anoms and not self.detrend:
            ts = self.ts_seas.loc[:,['anomalies']].copy(deep=True)
        if self.calc_anoms and self.detrend:
            ts = self.ts_seas.loc[:,['d_anomalies']].copy(deep=True)
        if not self.calc_anoms and self.detrend:
            ts = self.ts_seas.loc[:,['d_' + self.variable]].copy(deep=True)
        if not self.calc_anoms and not self.detrend:
            ts = self.ts_seas.loc[:,[self.variable]].copy(deep=True)

        labels=['WB','B','N','A','WA']

        sub, bins = pd.qcut(ts.iloc[:,0], 5, labels=labels, retbins=True)

        ts.loc[:,'cat'] = sub

        # if the flag qualitative is set to True (default is false)
        # then we search the years corresponding to the category
        if self.qualitative:
            self.value = str(self.value)
            if self.value not in labels:
                raise ValueError("category not in ['WB','B','N','A','WA']")
            else:
                tmp_df = ts[ts['cat'] == self.value].copy(deep=True)
                tmp_df.loc[:,'weights'] = 1. / len(tmp_df)
                self.analogs = tmp_df
                self.weights = self.analogs.loc[:,'weights'].values
                self.category = self.value
        # if value is quantitative we use the method ("quintiles" or "closest eight")
        else:
            self.value = float(self.value)
            if self.method == 'quintiles':
                bins[0] = -np.inf
                bins[-1] = np.inf
                category = labels[np.searchsorted(bins, np.float(self.value))-1]
                subset = ts[ts['cat'] == category]
                self.category = category
                tmp_df = subset.copy(deep=True)
                # calculates the weights (add to 1)
                self.analogs = self._calc_weights(tmp_df)
                self.quintiles = bins
            elif self.method == "closest 8":
                sub = (abs(self.value - ts.iloc[:,0])).sort_values()[:8].index
                tmp_df = ts.loc[sub,:].copy(deep=True)
                # calculates the weights (add to 1)
                self.analogs = self._calc_weights(tmp_df)
                self.category = self.analogs.loc[:,'cat'].values
            self.weights = self.analogs.loc[:,'weights'].values
        self.analog_years = self.analogs.index.year

    def proxy_repr(self, pprint=False, outfile=True):
        """
        proxy_dict is an OrderedDict
        """
        proxy_dict = od()
        proxy_dict['sitename'] = self.sitename
        proxy_dict['proxy_type'] = self.proxy_type
        proxy_dict['measurement'] = self.proxy_type
        proxy_dict['proxy_type'] = self.proxy_type
        proxy_dict['proxy_type'] = self.proxy_type
        proxy_dict['measurement'] = self.measurement
        proxy_dict['dating_convention'] = self.dating_convention
        proxy_dict['calendar'] = self.calendar
        proxy_dict['chronology'] = self.chronology

        proxy_dict['coords'] = self.coords
        proxy_dict['aspect'] = self.aspect
        proxy_dict['elevation'] = self.elevation

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
        if self.qualitative:
            proxy_dict['category'] = self.category
        else:
            if self.method == 'quintiles':
                proxy_dict['category'] = self.category
            elif self.method == 'closest 8':
                proxy_dict['category'] = ",".join(list(self.category))
        proxy_dict['analog_years'] = self.analog_years.tolist()
        proxy_dict['weights'] = list(self.weights)

        if pprint:
            pprint_od(proxy_dict)

        if outfile:
            # the name of the JSON file used to be created from the
            # proxy name
            # --------------------------------------------------------
            # proxy_name = self.sitename.replace(" ","_")
            # proxy_name = proxy_name.replace(".","")
            # fname = "{}.json".format(self.sitename.replace(" ","_"))
            # now the name of the JSON file is a parameter that is
            # passed to the script "proxy_oper" by the PHP layer
            # --------------------------------------------------------

            #with open(os.path.join(self.pjsons, fname),'w') as f:
            with open(os.path.join(self.pjsons, self.pfname),'w') as f:
                json.dump(proxy_dict, f)
        self.proxy_dict = proxy_dict

    def to_html(filename):
        if not(hasattr(self, 'analogs')):
            self.find_analogs()

    def plot(self):
        f = plot_season_ts(self)
        return f
