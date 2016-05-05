import os
import numpy as np
from numpy import ma
import json
import bottleneck  as bn
from matplotlib.mlab import detrend_linear
from scipy.stats import ttest_ind
from ..utils import seasons_params

try:
    import xarray as xray
except:
    try:
        import xray
    except ImportError:
        print('cannot import xarray or xray')

class analogs:
    """
    base class for analogs calculations,
    takes either a `proxy` or `ensemble` instance

    Parameters
    ----------

    dataset : string
            The name of the dataset for which the composite are to be calculated
            e.g. ['ersst','ncep','vcsn','gpcp']

    variable : string
            The name of the variable to consider for calculing the composite
            e.g. for 'ersst', it's 'sst', for 'ncep', can be ['hgt_1000',
            'hgt_850', etc]
    """
    def __init__(self, obj, dataset, variable):
        super(analogs, self).__init__()
        # the parent can be either an instance of a `proxy` or `ensemble` class
        self.parent = obj
        # the dataset to read from
        self.dataset = dataset
        # the variable to read
        self.variable = variable
        # the season is an attribute of the parent object
        # if the parent is an ensemble, the consistency has
        # already been checked
        self.season = self.parent.season
        self.analog_years = self.parent.analog_years
        # self.weights = self.parent.weights
        self.detrend = self.parent.detrend
        # wet get the location(s) (lon, lat) of the proxy
        # or ensemble of proxies here
        self.locations = {}
        if self.parent.description == 'proxy':
            self.locations[self.parent.sitename] = self.parent.extracted_coords
        elif self.parent.description == 'ensemble':
            for k in self.parent.dict_proxies.keys():
                self.locations[k] = self.parent.dict_proxies[k]['extracted_coords']

    def _read_dset_params(self):
        """
        reads in the Dataset parameters
        """
        with open(os.path.join(self.parent.djsons, 'datasets.json'), 'r') as f:
            dset_dict = json.loads(f.read())
        # dset_dict is a dictionnary holding useful metadata
        self.dset_dict = dset_dict[self.dataset][self.variable]

    def _check_domain(self, domain):
        if not(hasattr(self, 'dset_dict')):
            self._read_dset_params()
        self.domain = domain
        domain_dset = self.dset_dict['domain']
        if ( (self.domain[0] < domain_dset[0]) | (self.domain[1] > domain_dset[1])  \
            | (self.domain[2] < domain_dset[2]) | (self.domain[3] > domain_dset[3]) ):
            print("""ERROR! the domain for the composite is partly outside the limits of the dataset""")
            raise Exception("DOMAIN ERROR")

    def calculate_season(self):
        """
        calculates the season
        """
        self.seasons_params = seasons_params()

        if not(hasattr(self, 'dset_dict')):
            self._read_dset_params()

        # get the name of the file to open
        fname = self.dset_dict['path']

        # `dset` is now an attribute of the ensemble object
        self.dset = xray.open_dataset(fname)

        # get the variable and its index
        m_var = self.dset[self.variable].data
        index = self.dset['time'].to_index()

        # if the variable is rainfall, we calculate the running SUM
        if self.dset_dict['units'] in ['mm']:
            seas_field = bn.move_sum(m_var, self.seasons_params[self.season][0], \
                                          min_count=self.seasons_params[self.season][0], axis=0)
        # if not, then we calculate the running MEAN (average)
        else:
            seas_field = bn.move_mean(m_var, self.seasons_params[self.season][0], \
                                          min_count=self.seasons_params[self.season][0], axis=0)

        # get rid of the first nans in the time-series / fields after move_mean or move_sum
        seas_field = seas_field[(self.seasons_params[self.season][0]-1)::,:,:]
        index = index[(self.seasons_params[self.season][0]-1)::]

        # now selects the SEASON of interest
        iseas = np.where(index.month == self.seasons_params[self.season][1])[0]
        dates = index[iseas]
        seas_field = np.take(seas_field, iseas, axis=0)

        # if detrend is set to `True`, we detrend
        # detrend_linear from matplotlib.mlab is faster than detrend from scipy.signal
        if self.detrend:
            dseas_field = np.ones(seas_field.shape) * np.nan
            # if there is a mask, we have to test each variable
            if 'mask' in self.dset.data_vars:
                for ilat in range(dseas_field.shape[1]):
                    for ilon in range(dseas_field.shape[2]):
                        if np.logical_not(np.all(np.isnan(seas_field[:,ilat, ilon]))):
                            dseas_field[:,ilat, ilon] = detrend_linear(seas_field[:,ilat,ilon]) \
                            + seas_field[:,ilat,ilon].mean()

            # if not, we can proceed over the whole dataset
            else:
                for ilat in range(dseas_field.shape[1]):
                    for ilon in range(dseas_field.shape[2]):
                        dseas_field[:,ilat, ilon] = detrend_linear(seas_field[:,ilat,ilon]) \
                        + seas_field[:,ilat,ilon].mean()

            self.dset['dates'] = (('dates',), dates)
            self.dset['seas_var'] = (('dates', 'latitudes', 'longitudes'), dseas_field)

        # if detrend is False, then just add the seaosnal values
        else:
            self.dset['dates'] = (('dates',), dates)
            self.dset['seas_var'] = (('dates', 'latitudes', 'longitudes'), seas_field)

    def composite(self, climatology=(1981, 2010),  test=True, repeats=True, weighting=False):
        """
        calculates the composite anomalies (and the Student t-test)
        from the seasonal values

        Parameters
        ----------

        climatology : tuple (start_year, end_year)
                The climatological period

        test : Boolean (default = True)
                Whether to calculate the Student T-test (p-value)

        repeats : Boolean (default = False)
                whether to include the repeated years
                only applies when an `ensemble` object is passed to
                the analog class, where some years can be sampled
                repeatedly

        weigthing : Boolean (default = False)
                if True, calculate the composite anomaly
                weigthed by the inverse absolute difference between
                the proxy value and the analog seasons values
        """
        self.climatology = climatology

        # if we forgot to calculate the seasonal aggregate
        if not(hasattr(self, 'dset')):
            self.calculate_season()

        if repeats:
            # extract the composite sample: it INCLUDES the repeated years
            compos_s = xray.concat([self.dset['seas_var'].sel(dates=str(y)) for y in self.analog_years], dim='dates')
            ayears = self.analog_years
        else:
            # extract the composite sample EXCLUDING the repeated years
            compos_s = xray.concat([self.dset['seas_var'].sel(dates=str(y)) for y in np.unique(self.analog_years)], dim='dates')
            ayears = np.unique(self.analog_years)


        # calculating the climatology
        clim = self.dset['seas_var'].sel(dates=slice(str(self.climatology[0]), \
                                                     str(self.climatology[1])))

        # calculate the anomalies WRT the climatology
        compos_a = compos_s - clim.mean('dates')

        if weighting:
            # if weigthing, multiply by the weights (come from
            # either a proxy or an ensemble)

            # cast that into a xray.DataArray sharing the same axis than the composite anomalies
            weights_arr = xray.DataArray(np.array(self.parent.weights), dims=('dates'), coords={'dates':compos_a.dates.data})

            compos_a = (compos_a * weights_arr) / weights_arr.sum('dates')

        # get the composite anomalies into a DataArray
        compos_a_x = xray.DataArray(ma.masked_array(compos_a, np.isnan(compos_a)), dims=('years','latitudes','longitudes'),
                              coords={'years':ayears, 'latitudes':self.dset.latitudes, 'longitudes':self.dset.longitudes})
        # if test is True, then the standard Student t-test is calculated
        if test:
            t, pvalues = ttest_ind(compos_s.data, clim.data, axis=0)
            # pvalues contains the p-values, we can delete the Test statistics
            del(t)

        # we drop the time and dates dimensions, which has
        # for effect to drop all the variables that depend on them
        self.dset = self.dset.drop(('dates','time'))

        # store the anomalies and the composite anomalies
        # in the xray Dataset

        self.dset['years'] = (('years',), ayears)

        self.dset['composite_sample'] = compos_a_x

        self.dset['composite_anomalies'] = compos_a_x.mean('years')

        self.dset['weights'] = xray.DataArray(np.array(self.parent.weights), dims=('years'), coords={'years':compos_a_x.years.data})

        # saves the p-values
        self.dset['pvalues'] = \
        (('latitudes', 'longitudes'), pvalues)

        # set the attributes
        self.dset['latitudes'].attrs['units'] = 'degrees_north'
        self.dset['latitudes'].attrs['long_name'] = 'Latitudes'
        self.dset['latitudes'].attrs['axis'] = 'Y'
        self.dset['longitudes'].attrs['units'] = 'degrees_east'
        self.dset['longitudes'].attrs['long_name'] = 'Longitudes'
        self.dset['longitudes'].attrs['axis'] = 'X'
        self.dset['composite_sample'].attrs['missing_value'] = -999.9
        self.dset['composite_sample'].attrs['_FillValue'] = -999.9
        self.dset['composite_anomalies'].attrs['missing_value'] = -999.9
        self.dset['composite_anomalies'].attrs['_FillValue'] = -999.9

        return self

    def save_to_file(self, fname=None):
        nc = self.dset[['composite_sample','composite_anomalies', 'pvalues']]
        nc = nc.to_netcdf(fname)
        self.dset.close()

    def close(self):
        self.dset.close()
