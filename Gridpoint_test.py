import os
import re
from pyproj import Proj, transform
import time

from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')
from geopy.distance import distance
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import h5py
from scipy import stats

import shapefile
from pyhdf.SD import SD, SDC
from mpl_toolkits.basemap import Basemap
import pandas as pd
import xarray as xr

from matplotlib.colors import ListedColormap, BoundaryNorm
from datetime import datetime, timedelta

# 測試檔案
path = r"D:\Grid_poiint\GRIDCRO2D_d01_2023_03_02"
ncfile = nc.Dataset(path)
print(ncfile)

lon_EMIS = ncfile.variables['LON']
lat_EMIS = ncfile.variables['LAT']

print(lon_EMIS[0][0].shape)
print(lat_EMIS[0][0])



