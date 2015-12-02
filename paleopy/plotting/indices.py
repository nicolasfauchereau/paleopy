import os
import json
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import zscore, ttest_1samp
from ..utils import seasons_params

class indices():
    """
    defines an index
    """
    def __init__(self, obj, name=None):
        self.parent = obj
        self.analog_years = np.unique(obj.analog_years)
        self.season = obj.season
        self.name = name

    def _read_csv(self):
        s_pars = seasons_params()[self.season]
        with open(os.path.join(self.parent.djsons, 'indices.json'), 'r') as f:
            dict_json = json.loads(f.read())
        if self.name is not None:
            data = {}
            mat = pd.read_csv(dict_json[self.name]['path'], index_col=0, parse_dates=True)
            mat = pd.rolling_mean(mat, s_pars[0])
            mat.dropna(inplace=True)
            mat = mat[mat.index.month == s_pars[1]]
            mat = mat.apply(zscore)
            data[self.name] = mat
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

    def composite(self):
        if not(hasattr(self, 'data')):
            self._read_csv()
        if isinstance(self.data, pd.core.frame.DataFrame):
            compos = pd.concat([self.data.ix[str(y)] for y in self.analog_years])
        elif isinstance(self.data, dict):
            compos = {}
            for i, k in enumerate(self.data):
                df = self.data[k]
                compos[k] = pd.concat([df.ix[str(y)] for y in self.analog_years])
        self.compos = compos

    def _plot_df(self, df, ax, ax_n = 0, pval = None):
        """
        plots a boxplot for one Series or a DataFrame with
        one column (one index)
        """

        b = df.boxplot(ax = ax, widths=0.85, patch_artist=True)

        if df.mean().values < 0:
            plt.setp(b['boxes'],facecolor='b', edgecolor='b', alpha=0.5)
            plt.setp(b['fliers'], color='gray', marker='+',visible=True)
            plt.setp(b['medians'],color='k',linewidth=1)
            ax.set_ylim(-3, 3)
            ax.plot(1, df.mean().values, 's', color='b', ms=10)
        else:
            plt.setp(b['boxes'],facecolor='r', edgecolor='r', alpha=0.5)
            plt.setp(b['fliers'], color='gray', marker='+',visible=True)
            plt.setp(b['medians'],color='k',linewidth=1)
            ax.set_ylim(-3, 3)
            ax.plot(1, df.mean().values, 's', color='r', ms=10)
        ax.axhline(0, c='k')
        if ax_n == 0:
            ax.set_ylabel('std.', fontsize=14)
        [l.set_fontsize(13) for l in ax.xaxis.get_ticklabels()]
        [l.set_fontsize(13) for l in ax.yaxis.get_ticklabels()]
        if pval is not None:
            if pval <= 0.1:
                ax.set_title("p={:4.2f}\n$\mu$={:4.2f}$\sigma$={:4.2f}".format(pval, \
                df.mean().values[0], \
                df.std().values[0]), \
                fontsize=12)
            else:
                ax.set_title("p={:4.2f}\n$\mu$={:4.2f}\n$\sigma$={:4.2f}".format(pval, \
                df.mean().values[0], \
                df.std().values[0]), \
                fontsize=12)


    def plot(self):

        if not(hasattr(self, 'compos')):
            self.composite()

        if isinstance(self.compos, pd.core.frame.DataFrame):

            l = len(self.compos.columns)
            f, axes = plt.subplots(nrows=1, ncols=l, figsize=(l,6), sharey=True)
            f.subplots_adjust(wspace=0.0, left=0.15, bottom=0.05, top=0.87)
            axes = axes.flatten('F')
            for i, k in enumerate(self.compos.columns):

                df = self.compos[[k]]

                t, pval = ttest_1samp(df.values, 0)

                self._plot_df(df, axes[i], ax_n=i, pval=pval[0])

        elif isinstance(self.compos, dict):

            l = len(self.compos.keys())
            f, axes = plt.subplots(nrows=1, ncols=l, figsize=(l,6), sharey=True)
            f.subplots_adjust(wspace=0.0, left=0.15, bottom=0.05, top=0.87)
            axes = axes.flatten('F')
            for i, k in enumerate(self.compos.keys()):

                df = self.compos[k]

                t, pval = ttest_1samp(df.values, 0)

                self._plot_df(df, axes[i], ax_n=i, pval=pval[0])

        return f
