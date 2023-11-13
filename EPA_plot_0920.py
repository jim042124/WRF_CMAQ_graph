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

def EPA_PM(air_area, EPA_date_day): # 抓出EPA的數據(記得要設定excel的日期部分(儲存格))
    path = r"C:\Users\Jimmy\Desktop\Tropomi_data\No2_0703\空氣品質監測小時值資料(一般污染物,每日更新) (2023-03).csv"
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

    # 設定時間序列
    date_time = []
    start_date = datetime(2023, 3, EPA_date_day-1)
    end_date = datetime(2023, 3, EPA_date_day)
    
    current_date = start_date
    while current_date < end_date:
        current_date += timedelta(days=1)
        date_time.append(current_date.strftime('%Y-%m-%d'))
    print(date_time)
    data_pm = []
    # print(data_PM)

    for date in date_time:
        for siteid in data_station.loc[filt_area]['"siteid"']:
            filt = (data_PM['"siteid"'] == siteid) & (data_PM['"itemengname"'] == 'NO2') & (data_PM['"monitordate"'] == date)
            # 設定EPA測站的時間
            data_pm.append((data_PM[filt]['"monitorvalue13"']).values) 
            # print(data_pm.append((data_PM[filt]['"monitorvalue00"']).values))

    data_pm = [float(val) if val != 'x' else np.nan for val in data_pm]

    station_lonlat = [lon_station, lat_station, data_pm]
    # print(data_pm)
    return station_lonlat

monitordate = ["2023-03-01","2023-03-02","2023-03-03"]

hour_data = ['"monitorvalue00"','"monitorvalue01"','"monitorvalue02"','"monitorvalue03"',
              '"monitorvalue04"','"monitorvalue05"','"monitorvalue06"','"monitorvalue07"',
              '"monitorvalue08"','"monitorvalue09"','"monitorvalue10"','"monitorvalue11"',
              '"monitorvalue12"','"monitorvalue13"','"monitorvalue14"','"monitorvalue15"',
              '"monitorvalue16"','"monitorvalue17"','"monitorvalue18"','"monitorvalue19"',
              '"monitorvalue20"','"monitorvalue21"','"monitorvalue22"','"monitorvalue23"',]

path = r"C:\Users\Jimmy\Desktop\Tropomi_data\No2_0703\空氣品質監測小時值資料(一般污染物,每日更新) (2023-03).csv"
data_PM = pd.read_csv(path)
siteid = 3  #富貴角:84, 基隆:1, 萬里:3

def EPA_NO2(date, hour, species, siteid):
    filt = (data_PM['"siteid"'] == siteid) & (data_PM['"itemengname"'] == species) & (data_PM['"monitordate"'] == monitordate[date])
    return filt

NO2_data = []
NO_data = []
pm25_data = []

for date in range(len(monitordate)):
    for hour in range(len(hour_data)):
        filt = EPA_NO2(date, hour, 'NO2', 78)
        NO2_data.append(float((data_PM[filt][hour_data[hour]]).values[0])) # value[0]代表直接取出數值出來
        filt = EPA_NO2(date, hour, 'NO', 78)
        NO_data.append(float((data_PM[filt][hour_data[hour]]).values[0])) # value[0]代表直接取出數值出來
        filt = EPA_NO2(date, hour, 'PM2.5', 78)
        pm25_data.append(float((data_PM[filt][hour_data[hour]]).values[0])) # value[0]代表直接取出數值出來


hour_ticks = ["00","01","02","03","04","05","06","07","08","09",
              "10","11","12","13","14","15","16","17","18","19",
              "20","21","22","23"]

mouth_ticks = ["03/01","03/02","03/03"]

hour_title = []
for month_tick in mouth_ticks:
    for hour_tick in hour_ticks:
        hour_title.append(month_tick +" "+ hour_tick)



def NO2_plot():
    x_hour = np.arange(0,72,1)
    plt.rcParams['font.sans-serif'] = ['Microsoft JhengHei']   
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, ax1 = plt.subplots(figsize=(20, 4))
    # ax1.plot(x_hour, NO2_data, label="馬祖",color='darkcyan')
    ax1.plot(x_hour, NO_data, label="NO",color='darkgoldenrod')
    ax1.plot(x_hour, NO2_data, label="NO2",color='darkcyan')
    ax1.set_ylabel("濃度(μg/m3)")
    ax1.set_xlabel("時間")
    ax1.set_ylim(0, 15)
    plt.xticks(x_hour, hour_title, rotation=45,fontsize=6)
    plt.legend(loc='upper left')
    plt.grid(axis='y')

    ax2 = ax1.twinx()
    ax2.plot(x_hour, pm25_data, label="pM2.5",color='slateblue')
    ax2.set_ylim(0,80)
    ax2.set_ylabel("濃度(pM2.5)")
    plt.legend(loc='upper right')
    plt.savefig("C:\\Users\\Jimmy\\Desktop\\衛星排放源修正\\EPA_濃度繪圖\\"+"EPA_馬公",dpi=500)
    plt.show()

NO2_plot()


