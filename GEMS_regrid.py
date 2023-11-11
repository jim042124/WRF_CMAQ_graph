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
from scipy.interpolate import griddata

# 衛星的資料
path = r"D:\GEMS_data\2023_3\GK2_GEMS_L2_20230302_0545_NO2_FW_DPRO_ORI.nc"
ncfile = nc.Dataset(path)

print(ncfile.groups['Data Fields'].variables.keys())
lon = ncfile.groups['Geolocation Fields'].variables["Longitude"][0:1000] # 改緯度
lon_transpose = np.transpose(lon)[0:500] # 改經度
lon = np.transpose(lon_transpose)

lat = ncfile.groups['Geolocation Fields'].variables["Latitude"][0:1000]
lat_transpose = np.transpose(lat)[0:500]
lat = np.transpose(lat_transpose)

no2_data = ncfile.groups['Data Fields'].variables["ColumnAmountNO2Trop"][0:1000]
no2_data_transpose = np.transpose(no2_data)[0:500]
no2_data = np.transpose(no2_data_transpose)
no2_data[no2_data<0] =  np.nan

# 排放源的資料
path_EMIS = r"D:\2023Emission\TT_EM_2023060"
ncfile_EMIS = nc.Dataset(path_EMIS)
NO2_EMIS = ncfile_EMIS.variables["NO2"]

# 排放源的經度緯度
path_grid = r"D:\Grid_poiint\GRIDCRO2D_d01_2023_03_02"
ncfile = nc.Dataset(path_grid)
print(ncfile)
lon_EMIS = ncfile.variables['LON'][0][0]
lat_EMIS = ncfile.variables['LAT'][0][0]


print(NO2_EMIS[0][0].shape)

interpolated_no2 = griddata((lon.flatten(), lat.flatten()), no2_data.flatten(), 
                            (lon_EMIS, lat_EMIS), method='linear')

print(interpolated_no2.shape)

plt.pcolormesh(lon_EMIS, lat_EMIS, interpolated_no2)
plt.show()

# plt.pcolormesh(lon, lat, no2_data)
# plt.show()



# # 這邊是做內插的部分
# lon_target, lat_target = np.meshgrid(np.linspace(min(lon.flatten()), max(lon.flatten()), 200), 
#                                      np.linspace(min(lat.flatten()), max(lat.flatten()), 180))

# values_interp_no2 = griddata((lon.flatten(), lat.flatten()), no2_data.flatten(), (lon_target, lat_target), method='linear')
# print(values_interp_no2)

# print(lon_target.shape)
# print(lat_target.shape)
# print(values_interp_no2.shape)

taiwan_boundary = r'C:\Users\Jimmy\Desktop\python\WRF\taiwan_mpaboundary\COUNTY_MOI_1090820'
China_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\Chian_coastline\gadm36_CHN_1"

m = Basemap(projection='merc',
            llcrnrlat=20, urcrnrlat = 42, llcrnrlon=100, urcrnrlon = 130.25)   
#llcrnrlat=15, urcrnrlat = 40, llcrnrlon=90, urcrnrlon = 130.25
#llcrnrlat=21.58, urcrnrlat = 35.51, llcrnrlon=102.19, urcrnrlon = 124.73

m.readshapefile(taiwan_boundary,'COUNTY_MOI_1090820') 
m.readshapefile(China_coastline,'gadm36_CHN_1')                

# x, y = m(lons, lats)
# mask = (x > m.llcrnrlon) & (x < m.urcrnrlon) & (y > m.llcrnrlat) & (y < m.urcrnrlat)
# aod_masked = np.ma.masked_where(~mask, data_albedo_3)


# m.drawcoastlines()

lon_lines = [121]
lat_lines = [24]

# 繪製經緯度線
m.drawmeridians(lon_lines, labels=[False, False, False, True])
m.drawparallels(lat_lines, labels=[True, False, False, False])

# lon, lat = m(lon_aod, lat_aod)
# m.drawmapscale(lon=121.6, lat=22.038, lon0=122.25, lat0=21.8,
#             length=100, barstyle='fancy', fontsize=8)


# lon, lat = m(lon_target, lat_target)
# # m.quiver(lon, lat, values_interp_u, values_interp_v, scale=100, width=0.002)
# m.pcolormesh(lon, lat, values_interp_no2, cmap = 'jet')

lon, lat = m(lon_EMIS, lat_EMIS)
# m.quiver(lon, lat, values_interp_u, values_interp_v, scale=100, width=0.002)
m.pcolormesh(lon, lat, interpolated_no2*((100**2)/(6.02*(10**23)))*46*1000, cmap = 'jet',vmin=0, vmax=55)


# lon, lat = m(lon, lat)

# m.contourf(lon, lat, np.sqrt(wind_v**2 + wind_u**2), cmap='coolwarm')
# m.quiver(lon, lat, wind_u, wind_v, scale=100, width=0.002)
# m.pcolormesh(lon ,lat ,no2_data, cmap='jet')
plt.colorbar(label = 'NO2 (mg/m^2)',fraction = 0.1)

plt.title("GEMS_NO2_Regrid 0302")

plt.savefig("C:\\Users\\Jimmy\\Desktop\\衛星排放源修正\\"+"GEMS_0302_Regrid",dpi=500)
plt.show()



