import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
import shapefile
import sys
import os

from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
from pyproj import Proj, transformer
from scipy.interpolate import griddata
from matplotlib.font_manager import fontManager

import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['Microsoft JhengHei']
plt.rcParams['axes.unicode_minus'] = False

EPA_area = ["North","CM","Central","South","KP","Yilan","HD","TW"]
EPA_area = ["北部","竹苗","中部","南部","高屏","宜蘭","花東","臺灣"]
RMSE_fix = [20.22, 16.47, 15.17, 19.54, 18.60, 19.95, 12.98, 18.53]
RMSE_base = [21.89, 18.68, 17.49, 22.88, 19.61, 21.81, 14.05, 20.44]

bar_width = 0.35
x = np.arange(len(EPA_area))

plt.bar(x, RMSE_base, width=bar_width, label='RMSE_base', color='royalblue')
plt.bar(x + bar_width, RMSE_fix, width=bar_width, label='RMSE_fix', color='gold')

plt.xticks(x + bar_width / 2, EPA_area)
plt.title("RMSE_Model_Obs (20230302)")
plt.legend(loc = "lower right")
save_path = r"C:\\Users\\Jimmy\\Desktop\\衛星排放源修正\\"
plt.savefig(save_path + "RMSE_Model_Obs" ,dpi=500)
plt.show()
