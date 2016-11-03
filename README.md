# Paleopy

classes and methods to calculate and plot analog composites from
a paleoclimate proxy or an ensemble of paleoclimate proxies

- [Installing paleopy](#Installing-paleopy)
- [Calling paleopy from a script](#calling-paleopy-from-a-script)
- [Jupyter notebooks examples](#jupyter-notebooks-examples)

<hr>

## Installing paleopy

[paleopy](https://github.com/nicolasfauchereau/paleopy) is written in pure python, and can be installed by running

```
python setup.py install
```
in the top-level directory

## calling paleopy from a script  
<br>
The `scripts` folder contains two scripts illustrating respectively
how `paleopy` processes a *proxy* (`proxy_oper.py`) and an *ensemble
of proxies* (`ensemble_oper.py`) from user input passed to the command line.

+ **proxy_oper.py**

```
→ ./proxy_oper.py --help
usage: proxy_oper.py [-h] [-dj DJSONS] [-pj PJSONS] [-pn PFNAME] [-o OPATH]
                     [-n SITENAME] [-t PROXY_TYPE] [-lon LON] [-lat LAT]
                     [-dset DATASET] [-var VARIABLE] [-s SEASON] [-val VALUE]
                     [-q QUALITATIVE] [-per PERIOD] [-clim CLIMATOLOGY]
                     [-an CALC_ANOMS] [-dt DETREND] [-a ASPECT] [-e ELEVATION]
                     [-dc DATING_CONVENTION] [-cal CALENDAR] [-ch CHRONOLOGY]
                     [-m MEASUREMENT] [-v True|False]

optional arguments:
  -h, --help            show this help message and exit
  -dj DJSONS, --djsons DJSONS
                        the path to the jsons files defining the paths and
                        parameters of each dataset
  -pj PJSONS, --pjsons PJSONS
                        the path where to save the individual proxy json files
  -pn PFNAME, --pfname PFNAME
                        the name of the JSON file containing the information
                        for a single proxy
  -o OPATH, --opath OPATH
                        the path where to save the figures
  -n SITENAME, --name SITENAME
                        the name of the site
  -t PROXY_TYPE, --type PROXY_TYPE
                        the type of proxy (coral, Tree-ring, etc)
  -lon LON, --longitude LON
                        the longitude (decimal degree) of the proxy site
  -lat LAT, --latitude LAT
                        the latitude (decimal degree) of the proxy site
  -dset DATASET, --dataset DATASET
                        the dataset to interrogate to draw the analog years
  -var VARIABLE, --variable VARIABLE
                        the variable in the dataset to interrogate to draw the
                        analog years
  -s SEASON, --season SEASON
                        the season to which the proxy is sensitive
  -val VALUE, --value VALUE
                        the value for the proxy: can be either a float or a
                        string, if a string, must be in
                        ['WB','B','N','A','WA'] and the `qualitative` flag
                        must be set to True
  -q QUALITATIVE, --qualitative QUALITATIVE
                        a flag indicating whether the value passed (see above)
                        is qualitative or not, default to False: i.e.
                        interpret the value as a float
  -per PERIOD, --period PERIOD
                        the period from which to draw the analog seasons
  -clim CLIMATOLOGY, --climatology CLIMATOLOGY
                        the climatological period with respect to which the
                        anomalies are calculated
  -an CALC_ANOMS, --calc_anoms CALC_ANOMS
                        True if the anomalies are calculated, False otherwise.
                        Default is True
  -dt DETREND, --detrend DETREND
                        True if the time-series need detrended, False
                        otherwise. Default is True
  -a ASPECT, --aspect ASPECT
                        the aspect (in degrees, from 0 to 360)
  -e ELEVATION, --elevation ELEVATION
                        the elevation (in meters)
  -dc DATING_CONVENTION, --dating DATING_CONVENTION
                        the dating convention
  -cal CALENDAR, --calendar CALENDAR
                        the calendar year
  -ch CHRONOLOGY, --chronology CHRONOLOGY
                        the chronology control (i.e. 14C, Historic,
                        Dendrochronology, etc)
  -m MEASUREMENT, --measurement MEASUREMENT
                        the proxy measurement type (e.g. width for tree rings)
  -v True|False, --verbose True|False
                        Creates an output.json file with information about the
                        step and the percentage of progress. If ommited, default
                        values is False.
```

+ **ensemble_oper.py**

```
→ ./ensemble_oper.py --help
usage: ensemble_oper.py [-h] [-dj DJSONS] [-j PJSONS] [-o OPATH] [-s SEASON]
                        [-v True|False]

optional arguments:
  -h, --help            show this help message and exit
  -dj DJSONS, --djsons DJSONS
                        the path to the jsons files defining the paths and
                        parameters of each dataset
  -j PJSONS, --pjsons PJSONS
                        the directory containing the proxy json files
  -o OPATH, --opath OPATH
                        the directory in which to save the figures
  -s SEASON, --season SEASON
                        the season to consider: will be checked against the
                        individual proxies seasons
  -v True|False, --verbose True|False
                        Creates an output.json file with information about the
                        step and the percentage of progress. If ommited, default
                        values is False.
```

## Jupyter notebooks examples

In the `notebooks` folder, you will find 4 Jupyter notebooks:

+ [**`proxy.ipynb`**](https://github.com/nicolasfauchereau/paleopy/blob/master/notebooks/proxy.ipynb)
illustrates how a `proxy` (an individual *proxy*) class is instantiated and how
the methods are called to process it, including
the reconstruction of climate anomalies using the analog approach

+ [**`ensemble.ipynb`**](https://github.com/nicolasfauchereau/paleopy/blob/master/notebooks/ensemble.ipynb) illustrates how an `ensemble` (i.e. a collection of *proxies*) class is instantiated, then how ones reconstructs climate
anomalies using a network of proxies

+ [**`WR.ipynb`**](https://github.com/nicolasfauchereau/paleopy/blob/master/notebooks/WR.ipynb)
illustrates the reconstruction of Weather Regimes (WR) frequency anomalies from an instance of an `ensemble` class

+ [**`indices.ipynb`**](https://github.com/nicolasfauchereau/paleopy/blob/master/notebooks/indices.ipynb)
illustrates the reconstruction of anomalies for
a set of climate indices (currently the SOI, NINO 3.4 SSTs, the SAM index and the IOD index)

hi there 
