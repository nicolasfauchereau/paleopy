import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import addcyclic
import palettable

class scalar_plot:
    """
    scalar plot accepts a `analogs` object and implements
    methods to plot a scalar plot on a map
    """
    def __init__(self, analogs, center=None, robust=True, robust_percentile=1.0, \
              vmin=None, vmax=None, cmap=None, grid=True, domain=None, proj=None, res='c', test=None, border=True):

        self.analogs = analogs
        self.domain = domain
        self.center = center
        self.robust = robust
        self.robust_percentile = robust_percentile
        self.vmin = vmin
        self.vmax = vmax
        self.cmap = cmap
        self.grid = grid
        self.domain = domain
        self.proj = proj
        self.res = res
        self.test = test
        self.border = border


    def get_domain(self, domain):

        domain_dset = self.analogs.dset_dict['domain']

        if ( (self.domain[0] < domain_dset[0]) | (self.domain[1] > domain_dset[1])  \
            | (self.domain[2] < domain_dset[2]) | (self.domain[3] > domain_dset[3]) ):
            print("""ERROR! the domain for the analogsite is partly outside the limits of the dataset""")
            raise Exception("DOMAIN ERROR")

        else:
            # checks whether the first latitude is the northermost or southermost latitude
            latstart = self.analogs.dset['latitudes'].data[0]
            latend  = self.analogs.dset['latitudes'].data[-1]
            # if first latitude northermost, reverse the order of the domain selection
            if latstart > latend:
                dset_domain = self.analogs.dset.sel(latitudes=slice(self.domain[3], self.domain[2]), \
                                                   longitudes=slice(self.domain[0], self.domain[1]))
            # if not, go ahead with the domain limits as they are provided
            else:
                dset_domain = self.analogs.dset.sel(latitudes=slice(self.domain[2], self.domain[3]), \
                                                   longitudes=slice(self.domain[0], self.domain[1]))

            self.dset_domain = dset_domain

        return self


    def plot(self):

        """
        basemap plot of a scalar quantity
        """

        if self.domain:
            self.get_domain(self.domain)
        else:
            self.dset_domain = self.analogs.dset

        # get the lat and lons
        latitudes = self.dset_domain['latitudes'].data
        longitudes = self.dset_domain['longitudes'].data

        # get the data (analogsite anomalies)
        mat = self.dset_domain['composite_anomalies'].data
        pvalues = self.dset_domain['pvalues'].data

        if not(self.proj):
            self.proj = 'cyl'
        if self.proj in ['merc', 'cyl']:
            m = bm(projection=self.proj,llcrnrlat=latitudes.min(),urcrnrlat=latitudes.max(),\
                llcrnrlon=longitudes.min(),urcrnrlon=longitudes.max(),\
                lat_ts=0, resolution=self.res)
        if self.proj == 'moll':
            m = bm(projection='moll',lon_0=180, resolution=res)
        if self.proj in ['npstere','spstere']:
            m = bm(projection=self.proj,boundinglat=0,lon_0=0, resolution=self.res, round=True)
            # we add cyclic longitude
            mat, lon = addcyclic(mat, longitudes)
            pvalues, lon = addcyclic(pvalues, longitudes)
            longitudes = lon

        # m = pickle.load( open( "basemap.pickle", "rb" ) )

        # meshgrid lon and lat for plotting with basemap
        lons, lats = np.meshgrid(longitudes, latitudes)

        # ravel and removes nans for calculation of intervals etc
        calc_data = self.analogs.dset['composite_anomalies'].data
        calc_data = np.ravel(calc_data[np.isfinite(calc_data)])

        # the following is borrowed from xray
        # see: plot.py in xray/xray/plot
        if self.vmin is None:
            self.vmin = np.percentile(calc_data, self.robust_percentile) if self.robust else calc_data.min()
        if self.vmax is None:
            self.vmax = np.percentile(calc_data, 100 - self.robust_percentile) if self.robust else calc_data.max()

        del(calc_data)

        # Simple heuristics for whether these data should  have a divergent map
        divergent = ((self.vmin < 0) and (self.vmax > 0)) or self.center is not None

        # Now set center to 0 so math below makes sense
        if self.center is None:
            self.center = 0

        # A divergent map should be symmetric around the center value
        if divergent:
            vlim = max(abs(self.vmin - self.center), abs(self.vmax - self.center))
            self.vmin, self.vmax = -vlim, vlim

        # Now add in the centering value and set the limits
        self.vmin += self.center
        self.vmax += self.center

        # Choose default colormaps if not provided
        if self.cmap is None:
            if divergent:
                # get the colormap defined in the dset_dict
                self.cmap = eval(self.analogs.dset_dict['plot']['cmap'])
            else:
                # get the default matplotlib colormap
                self.cmap = plt.get_cmap()

        # ===============================================================
        # plots

        r, c = mat.shape

        f, ax = plt.subplots(figsize=(8*(c/r), 8), dpi=200)

        # the basemap instance gets attached to the created axes
        m.ax = ax

        # pcolormesh
        im = m.pcolormesh(lons, lats, ma.masked_invalid(mat), \
                          vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, latlon=True)

        # if test is defined, one contours the p-values for that level
        if self.test:
            #P = ma.where(self.dset['pvalues'].data <0.1 , 1, np.nan)
            #m.contourf(lon, lat, P, colors='0.99', \
            #           hatches=['...'], latlon=True, zorder=2, alpha=0.2)
            m.contour(lons, lats, pvalues, [self.test], latlon=True, \
                      colors='#8C001A', linewidths=1.5)

        # sets the colorbar
        cb = m.colorbar(im)
        [l.set_fontsize(14) for l in cb.ax.yaxis.get_ticklabels()]
        cb.set_label(self.analogs.dset_dict['units'],fontsize=14)

        # draw the coastlines
        m.drawcoastlines()

        # if one is plotting SST data, fill the continents
        if self.analogs.variable in ['sst','SST']:
            m.fillcontinents('0.8', lake_color='0.8')

        # if grid, plots the lat  / lon lines on the map
        if self.grid:
            # if it's hires ('f' or 'h'), every 10 degrees
            if self.res in ['h','f']:
                m.drawmeridians(np.arange(0, 360, 5), labels=[0,0,0,1], fontsize=14)
                m.drawparallels(np.arange(-80, 80+5, 5), labels=[1,0,0,0], fontsize=14)
            # if not hires, plots lines every 40 degrees
            else:
                m.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], fontsize=14)
                m.drawparallels(np.arange(-80, 80+10, 40), labels=[1,0,0,0], fontsize=14)

        if not(self.border):
            ax.axis('off')

        # set the title from the description in the dataset + variable JSON entry
        ax.set_title(self.analogs.dset_dict['description'], fontsize=14)

        return f

        self.dset_domain.close()

        plt.close(f)
