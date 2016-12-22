import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import addcyclic, shiftgrid
import palettable
from scipy.stats import scoreatpercentile

class vector_plot:
    """
    plots vector plots (uwnd + vwnd) on top of geopotential anomalies (ht)

    Parameters
    ----------

    ucompos : object, composite anomalies and metadata for zonal wind component
            a composite object coming from e.g. `uwnd_1000 = analogs(ens, 'ncep', 'uwnd_1000').composite()`

    vcompos : object, composite anomalies and metadata for meridional wind component
            a composite object coming from e.g. `vwnd_1000 = analogs(ens, 'ncep', 'vwnd_1000').composite()`

    hgtcompos : object, composite anomalies and metadata for HGT
            a composite object coming from e.g. `hgt_1000 = analogs(ens, 'ncep', 'hgt_1000').composite()`

    """
    def __init__(self, ucompos, vcompos, hgtcompos):
        
        """
        these are the composite objects
        """
        self.hgtcompos = hgtcompos
        self.ucompos = ucompos
        self.vcompos = vcompos
        
        """
        these are the datasets containing the composite anomalies and tests
        """      
        self.windspeed_dset = (self.ucompos.dset_compos ** 2 + self.vcompos.dset_compos ** 2).apply(np.sqrt)
        self.uwnd_dset  = ucompos.dset_compos
        self.vwnd_dset  = vcompos.dset_compos
        self.hgt_dset  = hgtcompos.dset_compos

    def _get_levels(self, data):
        """
        data can be either the data array attached to:

        + self.analogs.dset_compos['composite_anomalies']

        or one of the data arrays attached to

        + self.analogs.dset_compos['composite_sample']
        """

        # ravel and removes nans for calculation of intervals etc
        calc_data = np.ravel(data[np.isfinite(data)])

        # the following is borrowed from xray
        # see: plot.py in xray/xray/plot
        vmin = np.percentile(calc_data, 1)
        vmax = np.percentile(calc_data, 99)

        del(calc_data)

        # Simple heuristics for whether these data should  have a divergent map
        divergent = ((vmin < 0) and (vmax > 0))

        # A divergent map should be symmetric around the center value
        if divergent:
            center = 0
            vlim = max(abs(vmin), abs(vmax))
            vmin, vmax = -vlim, vlim

        if (vmin.dtype == 'float64') and (vmax.dtype == 'float64'):
            vmin = float("{:6.2f}".format(vmin))
            vmax = float("{:6.2f}".format(vmax))

        levels = np.linspace(vmin, vmax, num=10, endpoint=True)

        if divergent:
            neglevels = levels[levels < 0 ]
            poslevels = levels[levels > 0 ]
            levels = [neglevels,poslevels]
            # get the colormap defined in the dset_dict for HGT
            cmap = eval(self.hgtcompos.dset_dict['plot']['cmap'])
        else:
            # get the default matplotlib colormap
            cmap = plt.get_cmap()

        return vmin, vmax, levels, cmap
    
    def plot(self, domain = None, res='c', stepp=1, scale=8, test=0.1):
        
        """
        if the domain is actually defined, we select in lat and lon, making sure the 
        latitudes are increasing (from South to North)
        """
        
        if domain is not None: 
        
            latrev = (self.windspeed.latitudes[-1] < self.windspeed.latitudes[0])

            if latrev: 
                ugrid = self.uwnd_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
                vgrid = self.vwnd_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
                hgrid = self.hgt_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
                wgrid = self.windspeed_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
            else: 
                ugrid = self.uwnd_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
                vgrid = self.vwnd_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
                hgrid = self.hgt_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
                wgrid = self.windspeed_dset.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
        
        else: 
            
            ugrid = self.uwnd_dset
            vgrid = self.vwnd_dset
            hgrid = self.hgt_dset
            wgrid = self.windspeed_dset
            
        
        latitudes = wgrid.latitudes.data
        longitudes = wgrid.longitudes.data

        m = bm(projection='cyl',llcrnrlat=latitudes.min(),urcrnrlat=latitudes.max(),\
        llcrnrlon=longitudes.min(),urcrnrlon=longitudes.max(),\
        lat_ts=0, resolution=res, area_thresh=10000)

        lons, lats = np.meshgrid(longitudes, latitudes)
        
        # this colormap is quite good for wind speed
        # cmap = palettable.colorbrewer.sequential.Oranges_9.mpl_colormap
        

        """
        get the width and height of the figure
        """
        
        w = 10 
        
        h = np.ceil( (wgrid['composite_anomalies'].shape[0] / wgrid['composite_anomalies'].shape[1]) *  10 ) 
        
        f, ax = plt.subplots(figsize=(w,h))

        m.ax = ax
        
        m.drawmeridians(np.arange(0., 360. + 60, 60.), labels=[0,0,0,1], color='0.4', linewidth=0.5)
        m.drawparallels(np.arange(-80., 80. + 40., 40.), labels=[1,0,0,0], color='0.4', linewidth=0.5)
        
        """
        get the min, max, levels (for contours) and colormap
        """

        vmin, vmax, levels, cmap = self._get_levels(hgrid['composite_anomalies'].data)

        """
        plot using pcolormesh 
        """
        
        im = m.pcolormesh(lons, lats, hgrid['composite_anomalies'].data, cmap=cmap, vmin=vmin, vmax=vmax)
        
        cb = m.colorbar(im)

        cb.set_label("{}, {}".format(self.hgtcompos.dset_dict['short_description'], self.hgtcompos.dset_dict['units']), fontsize=14)

        ax.set_title("geopotential at 850 hPa \n{} and {}".format(self.ucompos.dset_dict['description'], \
                                                     self.vcompos.dset_dict['description']), fontsize=14)
        
        """
        plots the contours for HGT: not needed anymore
        """
        
#         if len(levels) == 2:
#             cn = m.contour(lons, lats, hgrid.data, levels = levels[0], cmap=plt.get_cmap('Blues'), linestyles='solid', latlon=True)
#             cp = m.contour(lons, lats, hgrid.data, levels = levels[1], cmap=plt.get_cmap('Reds'), latlon=True)
#         else: 
#             if hgrid.data.min() < 0: 
#                 cn = m.contour(lons, lats, hgrid.data, levels = levels[0], cmap=plt.get_cmap('Blues'), linestyles='solid', latlon=True)
#             elif hgrid.data.min() > 0: 
#                 cp = m.contour(lons, lats, hgrid.data, levels = levels[1], cmap=plt.get_cmap('Reds'), latlon=True)
                
        # plt.clabel(cn, fmt = '%4.0f', fontsize = 12, alpha=0.8, colors='k')
        # plt.clabel(cp, fmt = '%4.0f', fontsize = 12, alpha=0.8, colors='k')    
        
        """
        get the steps and plots the wind vectors ... stepp cannot really be determined 
        automatically and is therefore a parameter of the method `plotmap` of the class `vector_plot`
        """
        
        yy = np.arange(0, lats.shape[0], stepp)
        xx = np.arange(0, lons.shape[1], stepp)

        points = np.meshgrid(yy, xx)
        
        cmap_wind = palettable.colorbrewer.sequential.YlOrBr_9.mpl_colormap
            
#         Q = m.quiver(lons[points], lats[points], ugrid.data[points], vgrid.data[points], wgrid.data[points], \
#                      pivot='middle', scale=scale, cmap=plt.get_cmap('gray_r'), latlon=True)
    
        Q = m.quiver(lons[points], lats[points], ugrid['composite_anomalies'].data[points], vgrid['composite_anomalies'].data[points], \
                     pivot='middle', scale=scale, color='0.4', latlon=True)
        
        """
        test
        """
        
        mask = np.logical_or((c.ucompos.dset_compos['pvalues'].data < 0.1),(c.vcompos.dset_compos['pvalues'].data < test))
        
        ugrid_test = ma.masked_array(ugrid['composite_anomalies'].data, -mask) 
        vgrid_test = ma.masked_array(vgrid['composite_anomalies'].data, -mask) 
        
        QT = m.quiver(lons[points], lats[points], ugrid_test[points], vgrid_test[points], \
                     pivot='middle', scale=scale, color='k', latlon=True)        
        
        
        
#         m.streamplot(lons, lats, ugrid.data, vgrid.data, color='k', latlon=True, density=2.5, cmap=plt.cm.gray_r, linewidth=3*wgrid.data)
        
        l,b,w,h = ax.get_position().bounds
        
        """
        determine the wind vector key length as the 95th percentile of the wind speed 
        """

        wsq = np.round(scoreatpercentile(wgrid['composite_anomalies'].data, 95), decimals=1)
        
        qk = plt.quiverkey(Q, l+w-0.1, b-0.01, wsq, "{:4.2f} m/s".format(wsq), labelpos='E', fontproperties={'size':12}, coordinates='figure', zorder=20)

        m.drawcoastlines(color='0.2', linewidth=1)
        
        return f