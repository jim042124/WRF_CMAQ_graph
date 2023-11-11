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

def File_dir(path_dir): # 本程式碼主要運作於可將資料夾內的檔案全部讀取並使用
    file_paths = [] # 儲存的路徑
    for root, dirs, files in os.walk(path_dir): # os.walk返回出來的值有三個
        for file in files:
            file_path = os.path.join(root, file)
            file_paths.append(file_path)
    
    # 迭代並列印出檔案的路徑
    # for file_path in file_paths:
    #     print(file_path)
    return file_paths

def flatten(path): # 將衛星資料一維化並減少網格計算資源
    ncfile = nc.Dataset(path)
    # print(ncfile.groups['Data Fields'].variables.keys())
    lon = ncfile.groups['Geolocation Fields'].variables["Longitude"][:]
    lon_transpose = np.transpose(lon)[:]
    lon = np.transpose(lon_transpose)

    lat = ncfile.groups['Geolocation Fields'].variables["Latitude"][:]
    lat_transpose = np.transpose(lat)[:]
    lat = np.transpose(lat_transpose)

    no2_data = ncfile.groups['Data Fields'].variables["ColumnAmountNO2Trop"][:]
    no2_data_transpose = np.transpose(no2_data)[:]
    no2_data = np.transpose(no2_data_transpose)

    # 先設定需要的經度緯度
    lon_min, lon_max = 119.5, 122.25
    lat_min, lat_max = 21.8, 25.4

    # 取經度的索引(要經過轉換，矩陣先緯度後經度)
    lon_indices = []
    for i in range(len(lon.T)):
        select_data = lon.T[i][(lon.T[i] >= lon_min) & (lon.T[i] <= lon_max)]
        indices = np.where((lon.T[i] >= lon_min) & (lon.T[i] <= lon_max))
        if select_data != []:
            lon_indices.append(i)
    lon_indices = np.array(lon_indices)

    # 取緯度的索引    
    lat_indices = []
    for j in range(len(lat)):
        select_data = lat[j][(lat[j] >= lat_min) & (lat[j] <= lat_max)]
        indices = np.where(lat[j][(lat[j] >= lat_min) & (lat[j] <= lat_max)])
        if select_data != []:
            lat_indices.append(j)
    lat_indices = np.array(lat_indices)

    # 抓資料
    lon_new = []
    lat_new = []
    no2_new = []

    for x in lat_indices:
        lon_fix_lat = lon[x]
        lat_fix_lat = lat[x]
        no2_fix_lat = no2_data[x]
        for y in lon_indices:
            lon_fix_lon = lon_fix_lat[y]
            lat_fix_lon = lat_fix_lat[y]
            no2_fix_lon = no2_fix_lat[y]
            lon_new.append(lon_fix_lon)
            lat_new.append(lat_fix_lon)
            no2_new.append(no2_fix_lon)
    
    lon_new = np.array(lon_new)
    lat_new = np.array(lat_new)
    no2_new = np.array(no2_new)
    
    GEMS_flatten = [lon_new, lat_new, no2_new]
    return GEMS_flatten

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

def point_nearest(EPA_lonlat, GEMS_flatten): # 將衛星資料最靠近EPA的點取出
    diff_box = []
    air_mean = []
    no2_lon = []
    no2_lat = []
    no2_data = []

    for n in range(len(EPA_lonlat[0])):
        diff_box = []
        for i in range(len(GEMS_flatten[0])):
            diff = ((GEMS_flatten[0][i] - EPA_lonlat[0][n])**2 + (GEMS_flatten[1][i] - EPA_lonlat[1][n])**2)**(1/2)
            diff_box.append(diff)
        min_index = np.argmin(diff_box)
        no2_lon.append(GEMS_flatten[0][min_index])
        no2_lat.append(GEMS_flatten[1][min_index])
        no2_data.append(GEMS_flatten[2][min_index])
        print("目前已完成("+ str(n+1) +"/" + str(len(EPA_lonlat[0]))+")")
    
    save_Data = np.array([no2_lon, no2_lat, no2_data])
    return save_Data
    # np.savetxt('GEMS_output.txt', save_Data, delimiter=',')

def polyfit_graph(EPA_lonlat_2, no2_data): # 繪製趨勢線和相關性散佈圖的函式
    # 計算趨勢線及相關性 (正確的)
    m1, b1 = np.polyfit(EPA_lonlat_2, no2_data, 1)
    plt.plot(EPA_lonlat_2, m1*EPA_lonlat_2 + b1, color='green', label='Trend 1')

    corr = np.ma.corrcoef(EPA_lonlat_2, no2_data).data[0][1]
    print(corr)

    plt.scatter(EPA_lonlat_2, no2_data)
    plt.grid(linestyle=':')
    plt.xlabel('EPA_NO2',fontsize = 14)
    plt.ylabel('GEMS_NO2',fontsize = 14)

def graph_GEMS(path):
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
    no2_data = no2_data*((100**2)/(6.02*(10**23)))*46*1000
    print(ncfile.groups['Data Fields'].variables["ColumnAmountNO2Trop"])

    taiwan_boundary = r'C:\Users\Jimmy\Desktop\python\WRF\taiwan_mpaboundary\COUNTY_MOI_1090820'
    China_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\Chian_coastline\gadm36_CHN_1"
    Korea_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\korea_coastline\KOR_adm1"
    North_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\North_korea\PRK_adm0"

    m = Basemap(projection='merc',
                llcrnrlat=20, urcrnrlat = 42, llcrnrlon=100, urcrnrlon = 130.25)   
    # llcrnrlat=15, urcrnrlat = 40, llcrnrlon=90, urcrnrlon = 130.25
    # llcrnrlat=21.8, urcrnrlat = 25.4, llcrnrlon=119.5, urcrnrlon = 122.25
    # llcrnrlat=21.58, urcrnrlat = 35.51, llcrnrlon=102.19, urcrnrlon = 124.73
    # 東亞範圍：llcrnrlat=20, urcrnrlat = 42, llcrnrlon=100, urcrnrlon = 130.25
    # 論文範圍：llcrnrlat=21.51, urcrnrlat=25.51, llcrnrlon=119, urcrnrlon = 122.5
    m.readshapefile(taiwan_boundary,'COUNTY_MOI_1090820')
    m.readshapefile(China_coastline,'gadm36_CHN_1')              
    m.readshapefile(Korea_coastline,'KOR_adm1')
    m.readshapefile(North_coastline,'PRK_adm0')

    lon_lines = [121]
    lat_lines = [24]

    # 繪製經緯度線
    m.drawmeridians(lon_lines, labels=[False, False, False, True])
    m.drawparallels(lat_lines, labels=[True, False, False, False])


    lon, lat = m(lon, lat)
    # m.pcolormesh(lon, lat, no2_data, cmap = 'jet', vmin = 0, vmax = 20)

    m.pcolormesh(lon, lat, no2_data, cmap='jet',vmin= 0, vmax=70)
    # m.contourf(lon, lat, no2_data, cmap='jet')
    # m.contourf(lon, lat, no2_data, cmap='jet',levels=np.linspace(0, 55, 20))
    m.colorbar(label = 'UNIT:molec/cm2',fraction=0.09)
    plt.title("GEMS NO2_VCD 0302")

    #plt.savefig("C:\\Users\\Jimmy\\Desktop\\衛星排放源修正\\"+"GEMS_data_new_20230302",dpi=500)
    
def wind_field(ERA5_path):  # 用來繪製ERA5的檔案
    ncfile = nc.Dataset(ERA5_path)
    print(ncfile)
    print(ncfile.variables.keys()) # 查看檔案中的變數

    lon_wind, lat_wind = ncfile.variables['longitude'][:], ncfile.variables['latitude'][:]
    # print(lon_wind)
    # print(lat_wind)
    print(ncfile.variables['level'][:])

    lon, lat = np.meshgrid(lon_wind, lat_wind)
    wind_u, wind_v = ncfile.variables['u'][0][2], ncfile.variables['v'][0][2]

    taiwan_boundary = r'C:\Users\Jimmy\Desktop\python\WRF\taiwan_mpaboundary\COUNTY_MOI_1090820'
    China_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\Chian_coastline\gadm36_CHN_1"

    m = Basemap(projection='merc',
                llcrnrlat=21.51, urcrnrlat=25.51, llcrnrlon=119, urcrnrlon = 122.5)   
    # 論文中採用的經度緯度
    # llcrnrlat=21.51, urcrnrlat=25.51, llcrnrlon=119, urcrnrlon = 122.5
    # llcrnrlat=21.58, urcrnrlat = 35.51, llcrnrlon=102.19, urcrnrlon = 124.73

    m.readshapefile(taiwan_boundary,'COUNTY_MOI_1090820') 
    m.readshapefile(China_coastline,'gadm36_CHN_1')                

    lon_lines = [121]
    lat_lines = [24]

    # 繪製經緯度線
    m.drawmeridians(lon_lines, labels=[False, False, False, True])
    m.drawparallels(lat_lines, labels=[True, False, False, False])

    lon, lat = m(lon, lat)

    # m.contourf(lon, lat, np.sqrt(wind_v**2 + wind_u**2), cmap='jet', vmin = 0, vmax = 12)
    # cbar = plt.colorbar(label = 'Wind_unit (m/s)')
    # cbar.set_ticks(np.linspace(0, 12, 10))
    # cbar.set_ticklabels(np.arange(10))

    m.quiver(lon, lat, wind_u, wind_v, scale=100, width=0.003,color = 'white')

    plt.title("ERA_5 wind_field 0625")

    # plt.savefig("C:\\Users\\Jimmy\\Desktop\\pm_aod\\aod_pm25_plot\\"+'Himawari_8_band_1',dpi=500)

def WRF_data_graph(): # 可以畫出wrf裡的變數
    wrfout_path = r"D:\WRFout\wrfout_d01_2023-03-02_06_00_00"
    ncfile = nc.Dataset(wrfout_path)
    # 要放進def的變數部分
    wrf_time = 0 # 產圖時間
    layer = 0 # 代表第一層
    ####################
    print(ncfile)
    print(ncfile.variables.keys())
    print(ncfile.variables["Times"][wrf_time])
    
    # 製作時間表示
    word_time = str("")
    for i in range(len(ncfile.variables["Times"][wrf_time])):
        # print(ncfile.variables["Times"][0][i])
        time_word = ncfile.variables["Times"][wrf_time][i]
        word_time += str(time_word)[2] 
    print(word_time)

    # 取出變數
    # 經度緯度可以直接放在 m(lon,lat)就好
    lat = ncfile.variables["XLAT"][wrf_time][:]
    lon = ncfile.variables["XLONG"][wrf_time][:]
    print(lat.shape)
    print(lon.shape)
    u_wind = ncfile.variables["U"][wrf_time][:][layer]
    v_wind = ncfile.variables["V"][wrf_time][:][layer]
    print(u_wind.shape)
    print(v_wind.shape)

    u_wind = np.transpose(u_wind)[:len(lat[0])]
    u_wind = u_wind.T
    v_wind = v_wind[:len(lat)]
    ###################################
    
    taiwan_boundary = r'C:\Users\Jimmy\Desktop\python\WRF\taiwan_mpaboundary\COUNTY_MOI_1090820'
    China_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\Chian_coastline\gadm36_CHN_1"
    Korea_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\korea_coastline\KOR_adm1"
    North_coastline = r"C:\Users\Jimmy\Desktop\Tropomi_data\North_korea\PRK_adm0"
    
    m = Basemap(projection='merc',
                llcrnrlat=20, urcrnrlat = 42, llcrnrlon=100, urcrnrlon = 130.25)   
    # 論文中採用的經度緯度
    # llcrnrlat=21.51, urcrnrlat=25.51, llcrnrlon=119, urcrnrlon = 122.5
    # llcrnrlat=21.58, urcrnrlat = 35.51, llcrnrlon=102.19, urcrnrlon = 124.73

    m.readshapefile(taiwan_boundary,'COUNTY_MOI_1090820') 
    m.readshapefile(China_coastline,'gadm36_CHN_1')  
    m.readshapefile(Korea_coastline,'KOR_adm1')
    m.readshapefile(North_coastline,'PRK_adm0')
                  
    lon_lines = [121]
    lat_lines = [24]

    # 繪製經緯度線
    m.drawmeridians(lon_lines, labels=[False, False, False, True])
    m.drawparallels(lat_lines, labels=[True, False, False, False])
    
    # 將經度緯度改成圖可用的經緯度
    lon_map, lat_map = m(lon, lat)

    # 風場塗色版
    wind_total = np.sqrt(u_wind**2 + v_wind**2)
    print(wind_total.shape)
    # m.contourf(lon_map, lat_map, wind_total, cmap='jet') 
    m.pcolormesh(lon_map, lat_map, wind_total, cmap='jet') 
    m.colorbar()

    # 調整風場大小
    wind_interval = 5
    lon_map = lon_map[::wind_interval, ::wind_interval]
    lat_map = lat_map[::wind_interval, ::wind_interval]
    u_wind = u_wind[::wind_interval, ::wind_interval]
    v_wind = v_wind[::wind_interval, ::wind_interval]
    
    # 目前先採用此風速
    # m.quiver(lon_map, lat_map, u_wind, v_wind, scale=400, width=0.002,color = 'black')
    m.barbs(lon_map, lat_map, u_wind, v_wind, barbcolor ='black', pivot = 'middle'
            ,sizes = {'emptybarb':0.01}, length = 6)  
    
    # emptybarb 是調整風速低於5時產出的小圈圈的半徑大小
    # 'width':5 , barb_increments = {'half':2,'full':4,'flag':20}

    plt.title("WRFout_0302_1400")
    # plt.savefig("C:\\Users\\Jimmy\\Desktop\\衛星排放源修正\\"+"Wrfwind_0302",dpi=500)

    # plt.show()

# 設計修正函式，包括開啟檔案跟修正部分
# def Emission_fix():
    # 開啟排放源修正前的檔案
    # 開啟排放源修正前的CMAQ的NO2氣柱密度
    # 開啟排放源修正的跑CMAQ的NO2氣柱密度檔案
    # 開啟GEMS檔案，可以寫個函式做簡化
    # 計算B
    # 計算新的排放量
    # 做一個新的檔案並將新的排放量放入




