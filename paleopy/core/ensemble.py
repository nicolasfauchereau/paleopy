# Python packages imports
import os
from glob import glob
import json

class ensemble:
    """
    base class for an ensemble, just gathers all the metadata and data contained
    in a set of JSON proxy files.

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

    season : string
             The season

    analog_years : list
                   list of analog years ()

    detrend : boolean
              whether the time-series have been detrended
              or not

    dict_proxies : dictionnary
                   A dictionnary containing the individual JSON
                   proxy files

    """
    def __init__(self, djsons='../jsons', pjsons='../jsons/proxies', season='DJF'):
        super(ensemble, self).__init__()
        # type
        self.description = 'ensemble'
        # `jsons` is the path to the individual proxies JSON files
        self.djsons = djsons
        self.pjsons = pjsons
        self.season = season

        lfiles = glob(os.path.join(self.pjsons, "*.json"))
        self.analog_years = []
        self.weights = []
        self.detrend = []
        seasons = []
        climatologies = []

        self.dict_proxies = {}
        for json_file in lfiles:
            with open(json_file, 'r') as f:
                d = json.loads(f.read())
                self.analog_years.extend(d['analog_years'])
                self.weights.extend(d['weights'])
                self.detrend.append(d['detrend'])
                seasons.append(d['season'])
                climatologies.append(d['climatology'])
                self.dict_proxies[d['sitename']] = d

        assert len(set(seasons)) == 1, "seasons in json files must be identical"

        if list(set(seasons))[0] != self.season:
            print("""ERROR! season in json files does not match the season passed
            to the ensemble object""")
            raise Exception("SEASON ERROR")

        if (sum(self.detrend) != len(self.detrend)):
            print("""ERROR! detrend set to True for some proxies but False for others""")
            raise Exception("DETREND ERROR")
        else:
            self.detrend = self.detrend[0]

        if len(set(tuple(x) for x in climatologies)) > 1:
            print("""ERROR! the climatologies are not identical between the proxies""")
            raise Exception("VLIMATOLOGY ERROR")
        else:
            self.climatology = climatologies[0]
