# Python packages imports
import os
from glob import glob
import json
from collections import OrderedDict as od
import pandas as pd

class Ensemble:
    """
    base class for an ensemble, gathers all the metadata and data contained
    in a set of JSON proxy files.

    if the proxies are `consistent` (i.e. same season, same method [`detrend` or
    not `detrend`] have been used to determine the analog years, a flag [`proxies_consistent`]
    is set to True, if not, it is set to False, and the composite anomalies are first calculated
    *independently* for each proxy, *then* averaged together)

    Parameters
    ----------

    jsons : string
            The path to the directory containing the json files
    season : string
            The season

    Attributes
    ----------

    djsons : string
            The path to the directory containing the datasets json files

    pjsons : string
            The path to the directory containing the individual proxy json files

    """
    def __init__(self, djsons='../jsons', pjsons='../jsons/proxies'):
        # set the description to be an ensemble
        self.description = 'ensemble'
        # `djsons` is the path to the individual proxies JSON files
        self.djsons = djsons
        # `pjson`
        self.pjsons = pjsons

        lfiles = glob(os.path.join(self.pjsons, "*.json"))
        self.analog_years = []
        self.weights = []
        self.detrend = []
        self.season = []
        self.climatology = []

        # creates an ordered dictionnary to store the individual proxies
        # dictionnaries coming from the JSON files
        self.dict_proxies = od()
        for json_file in lfiles:
            with open(json_file, 'r') as f:
                # d is the full dictionnary with all the metadata
                d = json.loads(f.read())
                self.analog_years.append(d['analog_years'])
                self.weights.append(d['weights'])
                self.detrend.append(d['detrend'])
                self.season.append(d['season'])
                self.climatology.append(tuple(d['climatology']))
                # dict_proxies is a dictionnary containing all the metadata for the proxy
                self.dict_proxies[d['sitename']] = d
                # df_proxies is the pandas.DataFrame version of the dictionnary, exportable in e.g. cvs or HTML
                self.df_proxies = pd.DataFrame(self.dict_proxies).T

        """
        test the consistency of the seasons, method [detrend or not] and climatologies
        between all selected proxies
        """
        if (len(set(self.season)) == 1) and (len(set(self.season)) == 1) and (len(set(self.climatology)) == 1):
            print("""
            seasons and method [`detrend` or not] are consistent among proxies,
            the composite anomalies for each proxy can be calculated
            for the ensemble as a whole using the extended list of analog years
            """)
            self.proxies_consistent = 1
            self.season = self.season[0]
            self.detrend = self.detrend[0]
            self.climatology = self.climatology[0]
        else:
            print("""
            seasons and method [`detrend` or not] are inconsistent among proxies,
            the composite anomalies for EACH proxy can only be calculated
            INDEPENDENTLY, THEN combined ...
            """)
            self.proxies_consistent = 0
            self.season = self.season
            self.detrend = self.detrend
            self.climatology = self.climatology

    def __repr__(self):
        """
        the internal representation of an Ensemble object when called interactively
        """
        return """
        Ensemble made of {0} proxies
        ------------------
        proxies sitenames:\n{1:<10}""".format(len(self.df_proxies), "\n".join(self.df_proxies.sort_index().index.tolist()))
