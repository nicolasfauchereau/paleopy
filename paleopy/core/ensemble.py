# Python packages imports
import os
from glob import glob
import json
from collections import OrderedDict as od

class ensemble:
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
        super(ensemble, self).__init__()
        # type
        self.description = 'ensemble'
        # `jsons` is the path to the individual proxies JSON files
        self.djsons = djsons
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
                self.analog_years.extend(d['analog_years'])
                self.weights.extend(d['weights'])
                self.detrend.append(d['detrend'])
                self.season.append(d['season'])
                self.climatology.append(d['climatology'])
                self.dict_proxies[d['sitename']] = d

        """
        test the consistency of the seasons, method [detrend or not] and climatologies
        between all selected proxies
        """
        if (len(set(seasons)) == 1) and (len(set(seasons)) == 1) and (len(set(climatologies)) == 1):
            print("""
            seasons and method [`detrend` or not] are consistent among proxies
            """)
            self.proxies_consistent = 1
            self.season = self.seasons[0]
            self.detrend = detrend[0]
            self.climatology = climatologies[0]
        else:
            self.proxies_consistent = 0
            self.season = self.seasons
            self.detrend = detrend
            self.climatology = climatologies[0]
