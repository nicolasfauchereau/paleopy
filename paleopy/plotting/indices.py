import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.stats import zscore, ttest_1samp
from ..utils import seasons_params

class indices():
    """
    defines an index
    """
    def __init__(self, obj, jsons = './jsons/', name=None):
        self
        self.analog_years = np.unique(obj.analog_years)
        self.season = obj.season
        self.jsons = jsons
        self.name = name

    def _read_csv(self):
        s_pars = seasons_params()[self.season]
        with open(os.path.join(self.jsons, 'indices.json'), 'r') as f:
            dict_json = json.loads(f.read())
        if self.name is not None:
            data = pd.read_csv(dict_json[self.name]['path'], index_col=0, parse_dates=True)
            data = pd.rolling_mean(data, s_pars[0])
            data.dropna(inplace=True)
            data = data[data.index.month == s_pars[1]]
            data = data.apply(zscore)
        else:
            data = {}
            for k in dict_json.keys():
                mat = pd.read_csv(dict_json[k]['path'], index_col=0, parse_dates=True)
                mat = pd.rolling_mean(mat, s_pars[0])
                mat.dropna(inplace=True)
                mat = mat[mat.index.month == s_pars[1]]
                mat = mat.apply(zscore)
                data[k] = mat
        self.data = data
        return self

    def _plot_df(self, df, compos, ax, ax_n = 0, pval = None):

        b = compos.boxplot(ax = ax, widths=0.85, patch_artist=True)

        if compos.mean().values < 0:
            plt.setp(b['boxes'],facecolor='b', edgecolor='b', alpha=0.5)
            plt.setp(b['fliers'], color='gray', marker='+',visible=True)
            plt.setp(b['medians'],color='k',linewidth=1)
            ax.set_ylim(-3, 3)
            ax.plot(1, compos.mean().values, 's', color='b', ms=10)
        else:
            plt.setp(b['boxes'],facecolor='r', edgecolor='r', alpha=0.5)
            plt.setp(b['fliers'], color='gray', marker='+',visible=True)
            plt.setp(b['medians'],color='k',linewidth=1)
            ax.set_ylim(-3, 3)
            ax.plot(1, compos.mean().values, 's', color='r', ms=10)
        ax.axhline(0, c='k')
        if ax_n == 0:
            ax.set_ylabel('std.', fontsize=14)
        [l.set_fontsize(13) for l in ax.xaxis.get_ticklabels()]
        [l.set_fontsize(13) for l in ax.yaxis.get_ticklabels()]
        if pval is not None:
            if pval <= 0.1:
                ax.set_title("p={:4.2f}\n".format(pval), fontsize=12, fontweight='bold')
            else:
                ax.set_title("p={:4.2f}\n".format(pval), fontsize=12)


    def plot(self):

        if not(hasattr(self, 'data')):
            self._read_csv()

        if isinstance(self.data, pd.core.frame.DataFrame):

            compos = pd.concat([self.data.ix[str(y)] for y in self.analog_years])

            t, pval = ttest_1samp(compos.values, 0)

            f, ax = plt.subplots(figsize=(1.5,6))

            f.subplots_adjust(left=0.5)

            self._plot_df(self.data, compos, ax, ax_n=0, pval=pval[0])

        elif isinstance(self.data, dict):

            l = len(self.data)
            f, axes = plt.subplots(nrows=1, ncols=l, figsize=(l,6), sharey=True)
            f.subplots_adjust(wspace=0.0, left=0.15)
            axes = axes.flatten('F')
            for i, k in enumerate(self.data):
                df = self.data[k]

                compos = pd.concat([df.ix[str(y)] for y in self.analog_years])

                t, pval = ttest_1samp(compos.values, 0)

                self._plot_df(df, compos, axes[i], ax_n=i, pval=pval[0])

        return f
