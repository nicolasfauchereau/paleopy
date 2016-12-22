import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import addcyclic, shiftgrid
import palettable

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
        
        self.hgtcompos = hgtcompos
        self.ucompos = ucompos
        self.vcompos = vcompos
        self.uanoms  = ucompos.dset_compos['composite_anomalies']
        self.vanoms  = vcompos.dset_compos['composite_anomalies']
        self.hanoms  = hgtcompos.dset_compos['composite_anomalies'] 
        self.windspeed = np.sqrt(self.uanoms**2 + self.vanoms**2)

    def __get_levels(self, data):
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
            cmap = eval(hgtcompos.dset_dict['plot']['cmap'])
        else:
            # get the default matplotlib colormap
            cmap = plt.get_cmap()

        return vmin, vmax, levels, cmap
    
    def plot(self, domain = [0., 360., -90., 90.], res='c', stepp=1, scale=8):
        
        latrev = (self.windspeed.latitudes[-1] < self.windspeed.latitudes[0])
        
        if latrev: 
            ugrid = self.uanoms.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
            vgrid = self.vanoms.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
            hgrid = self.hanoms.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
            wgrid = self.windspeed.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[3],domain[2]))
        else: 
            ugrid = self.uanoms.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
            vgrid = self.vanoms.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
            hgrid = self.hanoms.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))
            wgrid = self.windspeed.sel(longitudes=slice(domain[0],domain[1]), latitudes=slice(domain[2],domain[3]))           
                
        latitudes = wgrid.latitudes.data
        longitudes = wgrid.longitudes.data

        m = bm(projection='cyl',llcrnrlat=latitudes.min(),urcrnrlat=latitudes.max(),\
        llcrnrlon=longitudes.min(),urcrnrlon=longitudes.max(),\
        lat_ts=0, resolution=res, area_thresh=10000)

        lons, lats = np.meshgrid(longitudes, latitudes)
        
        # this colormap is quite good for wind speed
        # cmap = palettable.colorbrewer.sequential.Oranges_9.mpl_colormap
        
        # this one is the one used for HGT
#         cmap = palettable.colorbrewer.diverging.RdYlBu_11_r.mpl_colormap
        cmap = eval(hgtcompos.dset_dict['plot']['cmap'])
        
        
        """
        get the width and height of the figure
        """
        
        w = 10 
        
        h = np.ceil( (wgrid.shape[0] / wgrid.shape[1]) *  10 ) 
        
        f, ax = plt.subplots(figsize=(w,h))

        m.ax = ax
        
        """
        get the min, max, levels (for contours) and colormap
        """

        vmin, vmax, levels, cmap = self.__get_levels(hgrid.data)

        """
        plot using pcolormesh 
        """
        
        im = m.pcolormesh(lons, lats, hgrid.data, cmap=cmap)
        
        cb = m.colorbar(im)

        cb.set_label('HG (m/s)', fontsize=14)

        """
        plots the contours for HGT 
        """
        
        if len(levels) == 2:
            cn = m.contour(lons, lats, hgrid.data, levels = levels[0], cmap=plt.get_cmap('Blues'), linestyles='solid', latlon=True)
            cp = m.contour(lons, lats, hgrid.data, levels = levels[1], cmap=plt.get_cmap('Reds'), latlon=True)
        else: 
            if hgrid.data.min() < 0: 
                cn = m.contour(lons, lats, hgrid.data, levels = levels[0], cmap=plt.get_cmap('Blues'), linestyles='solid', latlon=True)
            elif hgrid.data.min() > 0: 
                cp = m.contour(lons, lats, hgrid.data, levels = levels[1], cmap=plt.get_cmap('Reds'), latlon=True)
                
        # plt.clabel(cn, fmt = '%4.0f', fontsize = 12, alpha=0.8, colors='k')
        # plt.clabel(cp, fmt = '%4.0f', fontsize = 12, alpha=0.8, colors='k')    
        
        """
        get the steps and plots the wind vectors ... stepp cannot really be determined 
        automatically and is therefore a parameter of the method `plotmap` of the class `vector_plot`
        """
        
        yy = np.arange(0, lats.shape[0], stepp)
        xx = np.arange(0, lons.shape[1], stepp)

        points = np.meshgrid(yy, xx)
            
        Q = m.quiver(lons[points], lats[points], ugrid.data[points], vgrid.data[points], wgrid.data[points], \
                     pivot='middle', scale=scale, cmap=plt.cm.gray_r, latlon=True)

        l,b,w,h = ax.get_position().bounds

        qk = plt.quiverkey(Q, l+w-0.1, b-0.01, 0.25, "0.25 m/s", labelpos='E', fontproperties={'size':12}, coordinates='figure', zorder=20)

        m.drawcoastlines(color='0.4')
        
        return f