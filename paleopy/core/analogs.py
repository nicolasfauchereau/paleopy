import os
import numpy as np
from numpy import ma
import json
import xray
import bottleneck  as bn
from matplotlib.mlab import detrend_linear
from scipy.stats import ttest_ind

class analogs:
    """
    base class for analogs calculations,
    takes either a `proxy` or `ensemble` instance
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
        self.detrend = self.parent.detrend

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
            self.read_dset_params()
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
        seasons_params['Year (Jan. - Dec.)'] = (12, 12)
        seasons_params['Hydro. year (Jul. - Jun.)'] = (12, 6)
        self.seasons_params = seasons_params

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


    def composite(self, climatology=(1981, 2010),  test=True):
        """
        calculates the composite anomalies (and the Student t-test)
        from the seasonal values
        """
        self.climatology = climatology

        # we forgot to calculate the seasonal aggregate
        if not(hasattr(self, 'dset')):
            self.calculate_season()

        lyears = self.analog_years
        index = self.dset.dates.to_index()

        aseas = []
        # getting all the indices corresponding to years in 'analog_years'
        # for the season of interest, repeats are allowed
        for y in self.analog_years:
            aseas.append( np.where( (index.year == y) )[0][0] )

        # taking the composite
        self.compos = np.take(self.dset['seas_var'].data, aseas, axis=0)

        # calculating the climatology
        clim = self.dset['seas_var'].sel(dates=slice(str(self.climatology[0]), \
                                                     str(self.climatology[1])))

        self.clim = clim.data

        composite = self.compos.mean(0) - self.clim.mean(0)

        composite = ma.masked_array(composite, np.isnan(composite))

        # if test is True, then the standard Student t-test is calculated
        if test:
            t, pvalues = ttest_ind(self.compos, self.clim, axis=0)
            # pvalues contains the p-values, we can delete the Test statistics
            del(t)

        # if there is a mask, we multiply the composite anomalies by it
#         if 'mask' in self.dset.data_vars:
#             mask = self.dset['mask'].data
#             # saves the composite anomalies
#             self.dset['composite_anomalies'] = \
#             (('latitudes', 'longitudes'), composite * mask)
#             # saves the p-values
#             self.dset['pvalues'] = \
#             (('latitudes', 'longitudes'), pvalues)
#         else:
        # saves the composite anomalies
        self.dset['composite_anomalies'] = \
        (('latitudes', 'longitudes'), composite)
        # saves the p-values
        self.dset['pvalues'] = \
        (('latitudes', 'longitudes'), pvalues)

        return self

    def save_to_file(self, fname=None):
        nc = dset[['composite_anomalies', 'pvalues']]
        nc = nc.to_netcdf(fname)

    def close(self):
        """
        implements a `close` method to close the open datasets
        and quit cleanly
        """
        self.dset.close()
