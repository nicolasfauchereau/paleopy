import os
import numpy as np
from numpy import ma
import json
from copy import deepcopy
from itertools import chain
import xray
import bottleneck  as bn
from scipy.stats import ttest_ind
from ..utils import seasons_params, detrend_array

class Analogs:
    """
    base class for analogs calculations,
    takes either a `proxy` or `ensemble` instance
    """
    def __init__(self, obj, dataset, variable):
        # the parent can be either an instance of a `proxy` or `ensemble` class
        self.parent = obj
        # the dataset to read from
        self.dataset = dataset
        # the variable to read
        self.variable = variable
        # the seasons is an attribute of the parent object
        # if the parent is an ensemble, the consistency has
        # already been checked, and is contained in the .proxies_consistent attribute
        # if the ensemble if NOT consistent, `season` will be a list of lists
        self.season = self.parent.season
        self.analog_years = self.parent.analog_years
        # self.weights = self.parent.weights
        self.detrend = self.parent.detrend
        # wet get the location(s) (lon, lat) of the proxy
        # or ensemble of proxies here, in case it / they need to be plotted against the
        # anomalies
        self.locations = {}
        if self.parent.description == 'proxy':
            self.locations[self.parent.sitename] = self.parent.extracted_coords
        elif self.parent.description == 'ensemble':
            for k in self.parent.dict_proxies.keys():
                self.locations[k] = self.parent.dict_proxies[k]['extracted_coords']

    # todo: use the @property decorator on that one, because it is just an attribute

    def __read_dset_params(self):
        """
        reads in the Dataset parameters
        """
        with open(os.path.join(self.parent.djsons, 'datasets.json'), 'r') as f:
            dset_dict = json.loads(f.read())
        # dset_dict is a dictionnary holding useful metadata
        self.dset_dict = dset_dict[self.dataset][self.variable]

    # @property
    # def dset_dict(self):
    #     """reads in the Dataset parameters."""
    #     with open(os.path.join(self.parent.djsons, 'datasets.json'), 'r') as f:
    #         _dset_dict = json.loads(f.read())
    #     return self._dset_dict

    def __calculate_season(self, season):
        """
        calculates the seasonal values for the dataset and the variable
        passed as arguments when instantiating the `Analogs` Class

        at this stage, no detrend is applied even if it is set to True, it
        will be done if needed when calling the `__process_one_proxy` method
        """
        seasons_parameters = seasons_params()

        # if not(hasattr(self, 'dset_dict')):
        #     self.__read_dset_params()

        # `dset` is now an attribute of the ensemble object
        self.dset = xray.open_dataset(self.dset_dict['path'])

        self.dset = self.dset[[self.variable]]

        # get the variable and its index
        m_var = self.dset[self.variable].data
        index = self.dset['time'].to_index()

        # if the variable is rainfall, we calculate the running SUM over the length of the
        # season defined in seasons_parameters
        if self.dset_dict['units'] in ['mm']:
            seas_field = bn.move_sum(m_var, seasons_parameters[season][0], \
                                          min_count=seasons_parameters[season][0], axis=0)

        # if not, then we calculate the running MEAN (average) over the length of the
        # season defined in seasons_parameters
        else:
            seas_field = bn.move_mean(m_var, seasons_parameters[season][0], \
                                          min_count=seasons_parameters[season][0], axis=0)

        # get rid of the first nans in the time-series / fields after move_mean or move_sum
        seas_field = seas_field[(seasons_parameters[season][0]-1)::,:,:]
        index = index[(seasons_parameters[season][0]-1)::]

        # now selects the SEASON of interest
        iseas = np.where(index.month == seasons_parameters[season][1])[0]
        dates = index[iseas]
        seas_field = np.take(seas_field, iseas, axis=0)

        # creates a xarray Dataset which will contain the seasonal values
        dset_compos_dict = {}
        dset_compos_dict['latitudes'] = (('latitudes'), self.dset.latitudes.data)
        dset_compos_dict['longitudes'] = (('longitudes'), self.dset.longitudes.data)
        dset_compos_dict['dates'] = (('dates'), dates)
        dset_compos_dict['seas_var'] = (('dates','latitudes','longitudes'), seas_field)
        dset_compos = xray.Dataset(dset_compos_dict)
        del(dset_compos_dict)

        return dset_compos


    def __process_one(self, years, season, climatology=(1981, 2010), detrend=True, test=True, weighting=True, weights=None, repeats=False):
        """
        calculate the composite values for the dataset and variable passed as
        arguments of 'Analogs' for ONE proxy or ONE ensemble constituted from
        consistent proxies ...
        """

        if not(hasattr(self, 'dset_compos_dict')):
            dset_compos = self.__calculate_season(season)

        """
        first step, detrend if set to True
        """


        # if detrend is set to `True`, we detrend
        # detrend_linear from matplotlib.mlab is faster than detrend from scipy.signal
        if detrend:

            # we pop out the seasonal field
            seas_field = dset_compos['seas_var'].data

            # creates a numpy array containing NaNs with the same shape
            # dseas_field = np.ones(seas_field.shape) * np.nan

            # would be better to use e.g. `np.apply_along_axis`

            dseas_field = np.apply_along_axis(detrend_array, 0, seas_field)

            # # if there is a mask, we have to test each variable
            # if 'mask' in self.dset.data_vars:
            #     for ilat in range(dseas_field.shape[1]):
            #         for ilon in range(dseas_field.shape[2]):
            #             if np.logical_not(np.all(np.isnan(seas_field[:,ilat, ilon]))):
            #                 dseas_field[:,ilat, ilon] = detrend_linear(seas_field[:,ilat,ilon]) \
            #                 + seas_field[:,ilat,ilon].mean()
            #
            # # if no mask, we can proceed over the whole dataset
            # else:
            #     for ilat in range(dseas_field.shape[1]):
            #         for ilon in range(dseas_field.shape[2]):
            #             dseas_field[:,ilat, ilon] = detrend_linear(seas_field[:,ilat,ilon]) \
            #             + seas_field[:,ilat,ilon].mean()
            """
            over_write the variable `seas_var` in `dset_compos_dict` with the detrended version
            """

            dset_compos['seas_var'] = (('dates','latitudes','longitudes'), dseas_field)

        if repeats:
            # extract the composite sample: it INCLUDES the repeated years
            compos_s = xray.concat([dset_compos['seas_var'].sel(dates=str(y)) for y in years], dim='dates')
            ayears = years
        else:
            # extract the composite sample EXCLUDING the repeated years
            compos_s = xray.concat([dset_compos['seas_var'].sel(dates=str(y)) for y in np.unique(years)], dim='dates')
            ayears = np.unique(years)


        # calculating the climatology
        clim = dset_compos['seas_var'].sel(dates=slice(str(climatology[0]), \
                                                     str(climatology[1])))

        # calculate the anomalies WRT the climatology: composite anomalies
        compos_a = compos_s - clim.mean('dates')

        if weighting and weights is not None:
            # if weigthing, multiply by the weights (come from
            # either a proxy or an ensemble)

            # cast that into a xray.DataArray sharing the same axis than the composite anomalies
            weights_arr = xray.DataArray(np.array(weights), dims=('dates'), coords={'dates':compos_a.dates.data})

            compos_a = (compos_a * weights_arr) / weights_arr.sum('dates')

        # get the composite anomalies into a DataArray
        compos_a_x = xray.DataArray(ma.masked_array(compos_a, np.isnan(compos_a)), dims=('years','latitudes','longitudes'),
                              coords={'years':years, 'latitudes':dset_compos.latitudes, 'longitudes':dset_compos.longitudes})

        if test:
            t, pvalues = ttest_ind(compos_s.data, clim.data, axis=0)
            del(t)

        """
        get all that into a dictionnary, which we will transform into a xarray Dataset
        """

        dset_compos_dict = {}
        dset_compos_dict['years'] = ayears
        dset_compos_dict['latitudes'] = (('latitudes'), dset_compos.latitudes.data)
        dset_compos_dict['longitudes'] = (('longitudes'), dset_compos.longitudes.data)

        dset_compos_dict['composite_anomalies'] = (('latitudes','longitudes'), compos_a_x.mean('years').data)

        dset_compos_dict['pvalues'] = (('latitudes', 'longitudes'), pvalues)

        dset_compos_dict['composite_sample'] = (('years','latitudes','longitudes'), compos_a_x.data)

        if weighting:
            dset_compos_dict['weights'] = (('years', weights_arr.data))

        dset_compos = xray.Dataset(dset_compos_dict)

        # set the attributes
        dset_compos['latitudes'].attrs['units'] = 'degrees_north'
        dset_compos['latitudes'].attrs['long_name'] = 'Latitudes'
        dset_compos['latitudes'].attrs['axis'] = 'Y'
        dset_compos['longitudes'].attrs['units'] = 'degrees_east'
        dset_compos['longitudes'].attrs['long_name'] = 'Longitudes'
        dset_compos['longitudes'].attrs['axis'] = 'X'
        dset_compos['composite_sample'].attrs['missing_value'] = -999.9
        dset_compos['composite_sample'].attrs['_FillValue'] = -999.9
        dset_compos['composite_anomalies'].attrs['missing_value'] = -999.9
        dset_compos['composite_anomalies'].attrs['_FillValue'] = -999.9

        """
        returns the composite dataset, will be set as an attribute of `Analogs` later
        """
        return dset_compos


    def composite(self, weighting=True):

        """
        1st case: we are processing a proxy
        ___process_one(self, years, season, climatology=(1981, 2010), detrend=True, test=True, weighting=True, weights=None, repeats=False)
        """

        if self.parent.description == 'proxy':
            dset_compos = self.__process_one(self.analog_years,\
            self.season,\
            detrend = self.parent.detrend,\
            climatology = self.parent.climatology,\
            weighting = weighting,\
            weights = self.parent.weights)

            self.dset_compos = dset_compos

        """
        2nd case: we are processing a CONSISTENT ensemble of proxies: i.e. the sensitive season
        is the same among proxies, and all proxies share the same detrend attribute and climatologies
        """

        if (self.parent.description == 'ensemble') and self.parent.proxies_consistent == 1:
            """
            first flatten the list of lists of years and the list of lists of weights
            using itertools.chain(*list)
            """
            analog_years = list(chain(*self.analog_years))
            weights = list(chain(*self.parent.weights))

            season = self.season
            detrend = self.detrend
            climatology = self.parent.climatology

            dset_compos = self.__process_one(analog_years,\
            season,\
            detrend = detrend,\
            climatology = climatology,\
            weighting = weighting,\
            weights = weights,
            repeats = True)

            self.dset_compos = dset_compos

        """
        3rd case: the ensemble of proxies is inconsistent: i.e. the sensitive season
        are NOT the same among proxies, and all proxies do NOT share the same detrend attribute and climatologies

        Note that the length of e.g. the list of analog years do NOT need to be the same
        """

        if (self.parent.description == 'ensemble') and (self.parent.proxies_consistent == 0):
            """
            first determine the length of the ensemble
            """

            ensemble_length = len(self.parent.df_proxies)

            if any( [ensemble_length != len(self.analog_years), ensemble_length != len(self.parent.weights), ensemble_length != len(self.parent.climatology), ensemble_length != len(self.detrend)] ):
                print("the length of the ensemble {} is inconsistent with the length of one of the lists (analog_years, detrend, climatology, weights)".format(ensemble_length))
            else:
                l_dset_compos = []
                for i in range(ensemble_length):
                    years = self.analog_years[i]
                    season = self.season[i]
                    detrend = self.detrend[i]
                    climatology = self.parent.climatology[i]
                    weights = self.parent.weights[i]
                    dset_compos = self.__process_one(years, season, climatology=climatology, detrend=detrend, weights=weights, weighting=weighting)
                    l_dset_compos.append(dset_compos)

                # concatenate the datasets along a 'proxy' dimension
                self.dset_compos = xray.concat(l_dset_compos, dim='proxy')


    def save_to_file(self, fname=None):
        nc = self.dset_compos[['composite_anomalies', 'pvalues']]
        nc = nc.to_netcdf(fname)
        nc.close()

    def close(self):
        self.dset.close()
