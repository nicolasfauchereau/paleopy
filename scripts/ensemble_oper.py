#!/Users/nicolasf/anaconda/bin/python
import os
import sys
import argparse
import json
import matplotlib

matplotlib.use('Agg')

sys.path.insert(0, '../')

from paleopy import analogs
from paleopy import ensemble
from paleopy import WR
from paleopy.plotting import scalar_plot
from paleopy.plotting import indices
from paleopy.utils import save_progress
from matplotlib import pyplot as plt

"""
Global function for saving progress
"""

"""
parse command line arguments
"""

parser = argparse.ArgumentParser()

parser.add_argument('-dj', '--djsons', dest='djsons', type=str, default='../jsons', \
                    help='the path to the jsons files defining the paths and parameters of each dataset')

parser.add_argument('-j', '--pjsons', dest='pjsons', type=str, default='../jsons/proxies', \
                    help='the directory containing the proxy json files')

parser.add_argument('-o', '--opath', dest='opath', type=str, default='./outputs', \
                    help='the directory in which to save the figures')

parser.add_argument('-s', '--season', dest='season', type=str, default='DJF', \
                    help='the season to consider: will be checked against the individual proxies seasons')

parser.add_argument('-v', '--verbose', dest='verbose', type=bool, default=False, \
                    help='Output progress')

"""
goes from argparse Namespace to a dictionnary or key / value arguments
"""

vargs = vars(parser.parse_args())

"""
pop `opath` (the path where the outputs are saved) out of the dictionnary
"""

opath = vargs.pop('opath')

"""
pop `verbose` out of the dictionnary
"""
verbose = vargs.pop('verbose')

"""
instantiates an `ensemble` class, pass the `vargs` dict of keyword arguments to the class
"""
ens = ensemble(**vargs)

if verbose:
    save_progress(opath, 'Process the ensemble', 0)

"""
Creates output file array
"""

images = []

"""
instantiate the analog classes with the proxy for each dataset + variable we
want to map
"""

# ==============================================================================
"""
HGT 1000
"""
if verbose:
    save_progress(opath, 'HGT 1000 NZ', 20)

hgt_1000 = analogs(ens, 'ncep', 'hgt_1000').composite()
f = scalar_plot(hgt_1000, test=0.1, domain=[165, 180, -50., -30], proj='cyl', res='h', vmin=-25.0, vmax=25.0).plot()
f.savefig(os.path.join(opath, 'hgt_1000_NZ.png'))
plt.close(f)
images.append(
    {'id': 'NCEP_NZ_map_ensemble', 'title': 'New Zealand ensemble map for z1000', 'filename': 'hgt_1000_NZ.png'})

if verbose:
    save_progress(opath, 'HGT 1000 SH', 23)

f = scalar_plot(hgt_1000, test=0.1, proj='spstere', res='i', vmin=-25.0, vmax=25.0).plot()
f.savefig(os.path.join(opath, 'hgt_1000_SH.png'))
images.append({'id': 'NCEP_SH_map_ensemble', 'title': 'Southern Hemisphere ensemble map for z1000',
               'filename': 'hgt_1000_SH.png'})
plt.close(f)

if verbose:
    save_progress(opath, 'HGT 1000 SP', 26)

f = scalar_plot(hgt_1000, test=0.1, domain=[135, 290, -50., 10], proj='cyl', res='i', vmin=-25.0, vmax=25.0).plot()
f.savefig(os.path.join(opath, 'hgt_1000_SP.png'))
images.append(
    {'id': 'NCEP_SP_map_ensemble', 'title': 'Southwest Pacific ensemble map for z1000', 'filename': 'hgt_1000_SP.png'})
plt.close(f)

# plots the vector wind on top of the HGT1000

uwnd_1000 = analogs(ens, 'ncep', 'uwnd_1000').composite()
vwnd_1000 = analogs(ens, 'ncep', 'vwnd_1000').composite()

f = vector_plot(uwnd_1000, vwnd_1000, hgt_1000, ....).plot()



# ==============================================================================
"""
SST
"""

if verbose:
    save_progress(opath, 'SST NZ', 30)

sst = analogs(ens, 'ersst', 'sst').composite()
f = scalar_plot(sst, test=0.1, domain=[165, 180, -50., -30], proj='cyl', res='h', vmin=-1.5, vmax=1.5).plot()
f.savefig(os.path.join(opath, 'sst_nz.png'))
images.append({'id': 'SST_NZ_map_ensemble', 'title': 'New Zealand ensemble map for SSTa', 'filename': 'sst_nz.png'})
plt.close(f)

if verbose:
    save_progress(opath, 'SST SH', 33)

f = scalar_plot(sst, test=0.1, proj='spstere', res='i', vmin=-1.5, vmax=1.5).plot()
f.savefig(os.path.join(opath, 'sst_sh.png'))
images.append(
    {'id': 'SST_SH_map_ensemble', 'title': 'Southern Hemisphere ensemble map for SSTa', 'filename': 'sst_sh.png'})
plt.close(f)

if verbose:
    save_progress(opath, 'SST SP', 36)

f = scalar_plot(sst, test=0.1, domain=[135, 290, -50., 10], proj='cyl', res='i', vmin=-1.5, vmax=1.5).plot()
f.savefig(os.path.join(opath, 'sst_sp.png'))
images.append(
    {'id': 'SST_SP_map_ensemble', 'title': 'Southwest Pacific ensemble map for SSTa', 'filename': 'sst_sp.png'})
plt.close(f)

# ==============================================================================
"""
UWND at 850 and 200 hPa
"""

if verbose:
    save_progress(opath, 'UWND at 200hpa', 40)

uwnd = analogs(ens, 'ncep', 'uwnd_200').composite()
f = scalar_plot(uwnd, test=0.05, proj='spstere').plot()
f.savefig(os.path.join(opath, 'uwnd_200.png'))
images.append({'id': 'NCEP_SH_map_ensemble_200', 'title': 'Zonal Wind at 200hPa', 'filename': 'uwnd_200.png'})
plt.close(f)

if verbose:
    save_progress(opath, 'UWND at 850hpa', 45)
uwnd = analogs(ens, 'ncep', 'uwnd_850').composite()
f = scalar_plot(uwnd, test=0.05, proj='spstere').plot()
f.savefig(os.path.join(opath, 'uwnd_850.png'))
images.append({'id': 'NCEP_SH_map_ensemble_850', 'title': 'Zonal Wind at 850hPa', 'filename': 'uwnd_850.png'})
plt.close(f)

# ==============================================================================
"""
CLIMATE INDICES
"""

if verbose:
    save_progress(opath, 'Climate Indices', 50)

f = indices(ens).plot()
f.savefig(os.path.join(opath, 'indices_ensemble.png'))
images.append({'id': 'indices_ensemble', 'title': 'Climate Indices', 'filename': 'indices_ensemble.png'})
plt.close(f)

# ==============================================================================
"""
VCSN
"""
if verbose:
    save_progress(opath, 'VCSN Rain', 70)

vcsn = analogs(ens, 'vcsn', 'Rain').composite()
f = scalar_plot(vcsn, test=0.1, proj='cyl', res='h').plot(subplots=False)
f.savefig(os.path.join(opath, 'VCSN_rain_ensemble.png'))
images.append(
    {'id': 'VCSN_Rain_composite_anomaly', 'title': 'VCSN precipitation anomaly', 'filename': 'VCSN_rain_ensemble.png'})
plt.close(f)

if verbose:
    save_progress(opath, 'VCSN Mean Temperature', 71)

vcsn = analogs(ens, 'vcsn', 'TMean').composite()
f = scalar_plot(vcsn, test=0.1, proj='cyl', res='h', vmin=-1.5, vmax=1.5).plot(subplots=False)
f.savefig(os.path.join(opath, 'VCSN_temp_ensemble.png'))
images.append(
    {'id': 'VCSN_Tmean_composite_anomaly', 'title': 'VCSN Tmean anomaly', 'filename': 'VCSN_temp_ensemble.png'})
plt.close(f)

if verbose:
    save_progress(opath, 'VCSN Vapour Pressure', 72)

vcsn = analogs(ens, 'vcsn', 'Vapour_Pressure').composite()
f = scalar_plot(vcsn, test=0.1, proj='cyl', res='h', vmin=-1.5, vmax=1.5).plot(subplots=False)
f.savefig(os.path.join(opath, 'VCSN_VP_ensemble.png'))
images.append(
    {'id': 'VCSN_VP_composite_anomaly', 'title': 'VCSN Vapour Pressure anomaly', 'filename': 'VCSN_VP_ensemble.png'})
plt.close(f)

# ==============================================================================
"""
WEATHER REGIMES
"""

if ens.proxies_consistent == 0 and (len(set(self.season)) != 1):
    print("cannot process the WR plots if the proxies's seasons are inconsistent")
    pass
else:
    if verbose:
        save_progress(opath, 'Weather Regimes - Barplot', 80)

    w = WR(ens, classification='New Zealand')
    w.probs_anomalies(kind='many')
    f = w.plot_bar()
    f.savefig(os.path.join(opath, 'WR_ensemble.png'))
    images.append({'id': 'Ensemble_bar_plot', 'title': 'NZ weather regimes (Kidson Types) for ensemble composite',
                   'filename': 'WR_ensemble.png'})
    plt.close(f)

    if verbose:
        save_progress(opath, 'Weather Regimes - Heatmap', 90)

    f = w.plot_heatmap()
    f.savefig(os.path.join(opath, 'WR_heatmap.png'))
    images.append({'id': 'Ensemble_heat_map', 'title': 'Significance of regime changes for ensemble members',
                   'filename': 'WR_heatmap.png'})
    plt.close(f)




if verbose:
    save_progress(opath, 'Complete', 100)
# Save images list to json file
with open(os.path.join(opath, 'images.json'), 'w') as f:
    json.dump(images, f)
