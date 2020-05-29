"""

Program to create spectrum of spatial pH data 


"""

import numpy as np
import matplotlib.pyplot as plt
import random
from netCDF4 import Dataset
import pandas as pd
#
# read and plot pH model output data, from a netcdf file
# 
# define the netcdf reading function
def import_data(file_name):
    data_netcdf = Dataset(file_name, mode = 'r')
    data = {}
    for vname in list(data_netcdf.variables):
        data[str(vname)] = data_netcdf.variables[vname][:]
    data_netcdf.close()
    return data
#
path_in = '/Users/MMStoll/Python/Data/Ocean569_Data/WorldpH_Data/Surface_pH_1770_2000.nc'
pH_data = import_data(path_in)

test1 = pd.DataFrame(pH_data)
# plt.plot(pH_data['Year'],pH_data['pH'],color='firebrick')
# plt.grid()
# plt.xlabel('Months since Jan. 1960 (through 1995)')
# plt.ylabel('CO2 Levels (ppm)')
# plt.title('Mauna Loa Time Series')
# # plt.savefig(path_out1)

