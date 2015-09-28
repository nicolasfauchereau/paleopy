import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import addcyclic
import palettable

class vector_plot:
    def __init__(self, ucompos, vcompos):
        self.ucompos = ucompos
        self.vcompos = vcompos

        self.uanoms  = self.ucompos.dset['composite_anomalies']
        self.vanoms  = self.vcompos.dset['composite_anomalies']

        self.windspeed = np.sqrt(np.power(self.uanoms, 2) + np.power(self.vanoms, 2))

    def plotmap(self, domain = [0., 360., -90., 90.], res='c', stepp=2, scale=20):

        latitudes = self.windspeed.latitudes.data
        longitudes = self.windspeed.longitudes.data

        m = bm(projection='cyl',llcrnrlat=latitudes.min(),urcrnrlat=latitudes.max(),\
        llcrnrlon=longitudes.min(),urcrnrlon=longitudes.max(),\
        lat_ts=0, resolution=res)

        lons, lats = np.meshgrid(longitudes, latitudes)

        cmap = palettable.colorbrewer.sequential.Oranges_9.mpl_colormap

        f, ax = plt.subplots(figsize=(10,6))

        m.ax = ax

        x, y = m(lons, lats)

        im = m.pcolormesh(lons, lats, self.windspeed.data, cmap=cmap)

        cb = m.colorbar(im)
        cb.set_label('wind speed (m/s)', fontsize=14)

        Q = m.quiver(x[::stepp,::stepp], y[::stepp,::stepp], \
                     self.uanoms.data[::stepp,::stepp], self.vanoms.data[::stepp,::stepp], \
                     pivot='middle', scale=scale)

        l,b,w,h = ax.get_position().bounds

        qk = plt.quiverkey(Q, l+w-0.1, b-0.03, 5, "5 m/s", labelpos='E', fontproperties={'size':14}, coordinates='figure')

        m.drawcoastlines()

        return f
