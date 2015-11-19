#!/home/nicolasf/epd/bin/python                                                                                  
# -*- coding: utf-8 -*-
# ==================================================================================
# code reproject_Kidson_Types.py
# Description: reproject the kidson types in the daily DETRENDED 
# geopotential from NCEP
# created on 2012-10-12
# Nicolas Fauchereau <Nicolas.Fauchereau@gmail.com>
# ==================================================================================

### ==============================================================================================================
# some definitions here

### ==============================================================================================================
import os, sys
import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import netcdftime
from scipy.stats import zscore
from scipy.io.matlab import loadmat
from time_utils import remove_leap

dpath = '/home/nicolasf/data/'

### ==============================================================================================================
### some definitions
### ==============================================================================================================
def euclid(v1,v2):
    from scipy.stats import zscore
    '''Squared Euclidean Distance between two scalars or equally matched vectors
    
       USAGE: d = euclid(v1,v2)'''
    v1 = zscore(v1.flatten())
    v2 = zscore(v2.flatten())
    d2= np.sqrt(np.sum((v1-v2)**2))                                                                                                                                                                               
    return d2

def find_cluster(list_clus, dict_clus, X):
    """
    must return a list with 
    indice (time)
    name of attributed cluster
    position in list (indice)
    distance 
    difference with closest distance (separation index)
    """
    data_list = []
    for t in xrange(X.shape[0]):
        dlist = []
        for clus_name in list_clus:
            dlist.append(euclid(dict_clus[clus_name], X[t,...]))
            ranks = np.argsort(dlist)
            index = np.argmin(dlist)
        data_list.append([t, list_clus[index], index, dlist])
    return data_list

### ==============================================================================================================
### for comparisons with the "real" Kidson types classification 
from scipy.io.matlab import loadmat

kid_path = os.path.join(dpath + "/KidsonTypes/")

fname = "clusters_daily.mat"

matfile = loadmat(kid_path + fname, struct_as_record=False)

clusters = matfile['clusters'][0,0]

tclus = clusters.time

a = np.where( tclus[:,0] >= 1972)[0][0]
z = np.where( tclus[:,0] <= 2010 )[0][-1] + 1 

tclus = tclus[a:z,...]

# i12 = np.where(tclus[:,-1] == 12)[0]
# 
# tclus = tclus[i12,...]

### ==============================================================================================================
### name of the regimes 
name = clusters.name
name = name[a:z,...]

### makes the names and types flat for lookup 
names = []
for nm in name:
    names.append(str(nm[0][0]))
names = np.array(names)
del(name)

### ==============================================================================================================
### that above for comparison with the recalculated Kidson's types '

x = loadmat(os.path.join(dpath, "KidsonTypes", "h1000_clus.mat"), struct_as_record=False)
x = x['h1000']
x = x[0,0]
x.data.shape

# restrict the dataset to 1972 - 2010

a = np.where( (x.time[:,0] >= 1972) )[0][0]
z = np.where( (x.time[:,0] <= 2010) )[0][-1] + 1

x.time = x.time[a:z,...]
x.data = x.data[a:z,...]

from matplotlib.mlab import detrend_linear

### ==============================================================================================================
### detrend the data itself ?

datad = np.empty(x.data.shape)

for i in xrange(x.data.shape[1]):
    datad[:,i] = detrend_linear(x.data[:,i])

x.data = datad

clus_eof_file = loadmat(os.path.join(dpath, "KidsonTypes", "clus_eof.mat"), struct_as_record=False)
clus_eof = clus_eof_file['clus_eof'][0,0]

### normalize 
za = x.data -  np.tile(clus_eof.mean.T, (x.data.shape[0],1))

### multiply by the EOFs to get the Principal components 
pc = np.dot(za,clus_eof.vect)

pc_mean = clus_eof_file['pc_mean']

### normalize by the mean of the original PCs 
pc = pc - pc_mean

# detrend the PRINCIPAL COMPONENTS 
pcd = np.empty_like(pc)
for i in xrange(pc.shape[1]):
    pcd[:,i] = detrend_linear(pc[:,i])

### standardize by row 

pc = zscore(pc, axis=1)

### corrdinates of the 12 types on the original scores 
# clus_eof.clusmean.shape
# (12, 5)


### ==============================================================================================================
### from James's code 
# clusname={'TSW','T','SW','NE','R','HW','HE','W','HNW','TNW','HSE','H'};                                                                                                                                       
# regimes={{'TSW','T','SW','TNW'},{'W','HNW','H'},{'NE','R','HW','HE','HSE'}};
# regname={'Trough','Zonal','Blocking'};

list_clus = ['TSW','T','SW','NE','R','HW','HE','W','HNW','TNW','HSE','H']

dict_clus = {}
for i, k in enumerate(list_clus):
    dict_clus[k] = clus_eof.clusmean[i,...]

### ==============================================================================================================
### find the closest cluster in EOF space 
data_list = find_cluster(list_clus, dict_clus, pc)
data_listd = find_cluster(list_clus, dict_clus, pcd)

### get the names only
cluster_names_recalc = [data_list[i][1] for i in xrange(data_list.__len__())]
cluster_names_recalcd = [data_listd[i][1] for i in xrange(data_listd.__len__())]

### and see if it matches to the ones calculated previously by James 
matches = []
for i in xrange(len(names)):
    if names[i] == cluster_names_recalcd[i]:
        matches.append(1)
    else:
        matches.append(0)

matches.count(1)
matches.count(0)

clim_kid_rec = [ np.float(cluster_names_recalcd.count(nm)) / cluster_names_recalcd.__len__() for nm in list_clus]
clim_kid_orig = [ np.float(names.tolist().count(nm)) / names.tolist().__len__() for nm in list_clus]

### ==============================================================================================================
### plot the CLIMATOLOGICAL distribution of kidson types for the given SEASON 
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title("Kidson types distribution")
# ax.bar(np.arange(0, 12), np.array(clim_kid_orig), color="0.8", width=1.)
# plt.axvline(4,color="r", linewidth=2)
# plt.axvline(7,color="r", linewidth=2)
# plt.text(0.08,0.9,"Trough",transform = ax.transAxes, bbox=dict(facecolor='w', alpha=0.5), fontsize=18)
# plt.text(0.38,0.9,"Zonal",transform = ax.transAxes, bbox=dict(facecolor='w', alpha=0.5), fontsize=18)
# plt.text(0.7,0.9,"Blocking",transform = ax.transAxes, bbox=dict(facecolor='w', alpha=0.5), fontsize=18)
# ax.set_xticks(np.arange(0.5, 12.5))
# ax.set_xticklabels(list_clus, rotation="vertical", size='small')
# ax.set_ylim(0, max(clim_kid_orig))
# plt.ylabel("%")
# plt.grid()
# plt.savefig(os.path.join(fpath,"Kidson_types_clim_distrib_"+season+".png"),dpi=300)
# plt.close()

### ==============================================================================================================
### save the clusters 

### select one value per day (12 UTC)

### ==============================================================================================================
### indice for 12:00 UCT
i12 = np.where(tclus[:,-1] == 12)[0]

### select only 12 UTC, so one regime per day !
tclus = tclus[i12,...]

cluster_names_recalcd = np.array(cluster_names_recalcd)

cluster_names_recalcd = cluster_names_recalcd[i12,]

calendar, ifeb = remove_leap(tclus)

cluster_names_recalcd = np.delete(cluster_names_recalcd, ifeb)

data = zip(tclus, cluster_names_recalcd)

# f = open("/home/nicolasf/research/NIWA/paleo/data/cluster_names_recalcd.txt", "w")
# 
# for l in data:
#     f.write( str(l[0]).strip('[]') + " " + l[1]  + "\n")
# f.close()
# 
# ess = np.loadtxt("/home/nicolasf/research/NIWA/paleo/data/cluster_names_recalcd.txt", dtype={'names': ('years', 'month', 'day', 'time', 'regime'),'formats': ('i4', 'i4', 'i4', 'i4', 'S4')})
