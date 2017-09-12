

# readme to create conda environment and install paleopy 

## create a conda environment with necessary packages 

**1) create the environment**, and install the packages needed

**NOTE:** `conda` should be in `/data/pict/anaconda/bin`

```
conda create -n paleopy python=3.6 numpy scipy pandas matplotlib basemap netcdf4 h5py bottleneck dask xarray pip
```

**2) install palettable **

activate the environment `paleopy`

```
source activate paleopy
```

then install palettable

```
pip install palettable
```

## Install paleopy

clone or pull the master branch of the NIWA repo 

```
git clone https://github.com/niwa/paleopy.git
```

cd in the `paleopy` directory (where the file `setup.py`) is present

**IMPORTANT:**

The line 5 in `core/analogs.py` (see https://github.com/niwa/paleopy/blob/master/paleopy/core/analogs.py)

needs to be changed from 

```python
import xray 
```

to 

```python
import xarray as xray
```

prior to installation !

When this line has been changed, you can install paleopy in the `paleopy` environment 

```
source activate paleopy
```

```
python setup.py install
```



