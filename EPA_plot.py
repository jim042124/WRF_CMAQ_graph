import os
import re
from pyproj import Proj, transform
import time
import sys

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

# EPA_time = ["2022-10-13","2022-10-14","2022-10-15","2022-10-16","2022-10-17","2022-10-18"]
EPA_time = ["2021-04-17","2021-04-18","2021-04-19","2021-04-20","2021-04-21","2021-04-22","2021-04-23","2021-04-24","2021-04-25"]
# EPA_time = ["2023-05-28","2023-05-29","2023-05-30","2023-05-31"]

Monitor_time = ['"monitorvalue00"','"monitorvalue01"','"monitorvalue02"','"monitorvalue03"',
            '"monitorvalue04"','"monitorvalue05"','"monitorvalue06"','"monitorvalue07"',
            '"monitorvalue08"','"monitorvalue09"','"monitorvalue10"','"monitorvalue11"',
            '"monitorvalue12"','"monitorvalue13"','"monitorvalue14"','"monitorvalue15"',
            '"monitorvalue16"','"monitorvalue17"','"monitorvalue18"','"monitorvalue19"',
            '"monitorvalue20"','"monitorvalue21"','"monitorvalue22"','"monitorvalue23"']

air_area = ['北部空品區','雲嘉南空品區','竹苗空品區','中部空品區','高屏空品區','花東空品區','宜蘭空品區']

def EPA_PM(air_area, EPA_date_day, Monitor_Hour): # 抓出EPA的數據(記得要設定excel的日期部分(儲存格))
    path = r"C:\Users\Jimmy\Desktop\數值模擬\空氣品質監測小時值資料(一般污染物,每日更新) (2021-04).csv"
    data_PM = pd.read_csv(path)

    path_station = r"C:\Users\Jimmy\Desktop\空氣品質監測站基本資料.csv"
    data_station = pd.read_csv(path_station)

    # filt = (data_PM['"sitename"'] == '中壢') & (data_PM['"itemengname"'] == 'PM2.5')
    filt_area = data_station['"areaname"'].isin(air_area)
    # ['北部空品區','雲嘉南空品區','竹苗空品區','中部空品區','高屏空品區','花東空品區','宜蘭空品區']
    lon_station = np.array(data_station.loc[filt_area]['"twd97lon"'])
    lat_station = np.array(data_station.loc[filt_area]['"twd97lat"'])
    # station_lonlat = [lon_station, lat_station]
    # print(station_lonlat)
    # print(len(station_lonlat[0]))
    
    data_pm = []
    data_wind_speed = []
    data_wind_direction = []
    # print(data_PM)

    for siteid in data_station.loc[filt_area]['"siteid"']:
        filt = (data_PM['"siteid"'] == siteid) & (data_PM['"itemengname"'] == 'PM2.5') & (data_PM['"monitordate"'] == EPA_date_day)
        data_pm.append((data_PM[filt][Monitor_Hour]).values)
        filt_wind_speed = (data_PM['"siteid"'] == siteid) & (data_PM['"itemengname"'] == 'WIND_SPEED') & (data_PM['"monitordate"'] == EPA_date_day)
        data_wind_speed.append((data_PM[filt_wind_speed][Monitor_Hour]).values)
        filt_wind_direction = (data_PM['"siteid"'] == siteid) & (data_PM['"itemengname"'] == 'WIND_DIREC') & (data_PM['"monitordate"'] == EPA_date_day)
        data_wind_direction.append((data_PM[filt_wind_direction][Monitor_Hour]).values)
    



    # for date in date_time:
    #     for siteid in data_station.loc[filt_area]['"siteid"']:
    #         filt = (data_PM['"siteid"'] == siteid) & (data_PM['"itemengname"'] == 'NO2') & (data_PM['"monitordate"'] == date)
    #         # 設定EPA測站的時間
    #         data_pm.append((data_PM[filt]['"monitorvalue13"']).values) 
    #         # print(data_pm.append((data_PM[filt]['"monitorvalue00"']).values))

    data_pm = [float(val) if val != 'x' else np.nan for val in data_pm]
    data_wind_speed = [float(val) if val != 'x' else np.nan for val in data_wind_speed]
    data_wind_direction = [float(val) if val != 'x' else np.nan for val in data_wind_direction]

    station_lonlat = [lon_station, lat_station, data_pm, data_wind_speed, data_wind_direction]
    # print(data_pm)
    return station_lonlat

def EPA_point(EPA_lon, EPA_lat, PM25, wind_speed, wind_direction):
    taiwan_boundary = r'C:\Users\Jimmy\Desktop\python\WRF\taiwan_mpaboundary\COUNTY_MOI_1090820'

    map = Basemap(projection='merc',
                  resolution='i', fix_aspect=True,
                  llcrnrlat=21.8, urcrnrlat = 25.4, llcrnrlon=119.5, urcrnrlon = 122.25)
    # 讀取鄉政資料(邊界)
    map.readshapefile(taiwan_boundary,'COUNTY_MOI_1090820')
    # 設置地圖屬性
    # map.drawcoastlines(linewidth=1) # 海岸線的寬度
    map.drawparallels(np.arange(23.5,25), labels=[1,0,0,0], fontsize=10)
    map.drawmeridians(np.arange(120,122), labels=[0,0,0,1], fontsize=10)

    
    # 在地圖上標示位置
    x,y = map(EPA_lon, EPA_lat)
    color = ('#FFFFFF','#BBFFBB','#93FF93','#79FF79','#FFF0AC','#FFE66F','#FFBB00','#FFA600','#EA7500','#FF6600','#FF4D00','#CE0000','#930000','#4D0000','#AE00AE','#6F00D2','#000093','#28004D')
    levels = [0, 2, 7, 10, 15, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 400]
    my_cmap = ListedColormap(color)
    norm = BoundaryNorm(levels, len(levels))
    map.scatter(x, y, c=PM25, cmap=my_cmap, norm = norm, marker='o', s =20)

    wind_speed = np.full((len(wind_direction)), 10)
    wind_direction = wind_direction + np.full((len(wind_direction)), 180)
    u = [ws * np.sin(np.radians(wd)) for ws, wd in zip(wind_speed, wind_direction)]
    v = [ws * np.cos(np.radians(wd)) for ws, wd in zip(wind_speed, wind_direction)]
    map.quiver(x,y,u,v)
    
    # x,y = map(air_lon, air_lat)
    # map.scatter(x, y, c=air_statePM, cmap=my_cmap, norm = norm, marker='o', s =7)


    plt.title("Taiwan_Map")
    plt.xlabel('lon',fontsize=12,x=1)
    plt.ylabel("lat",fontsize=12,y=1,rotation=0)
    plt.colorbar()
    # plt.title("Taiwan_PM2.5 & Wind Field")
    # plt.savefig('C:\\Users\\Jimmy\\Desktop\\airbox\\plc\\' + '0315_13am5', dpi=300)
    # plt.show()
    # plt.clf()

for EPA_date_day in EPA_time:
    for Monitor_Hour, Hour in zip(Monitor_time, np.arange(0,23,1)):
        EPA_lon, EPA_lat, PM25, wind_speed, wind_direction = EPA_PM(air_area, EPA_date_day, Monitor_Hour)
        EPA_point(EPA_lon, EPA_lat, PM25, wind_speed, wind_direction)
        plt.title(str(EPA_date_day)+"_"+str(Hour))
        plt.savefig('C:\\Users\\Jimmy\\Desktop\\數值模擬\\2021_04\\' + str(EPA_date_day) +"_"+ str(Hour), dpi=500)
        plt.clf()
        print("Successful_" + str(EPA_date_day) +"_"+ str(Hour))
        print(wind_direction)
        print(wind_speed)

        







