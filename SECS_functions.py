# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 20:18:36 2023

@author: may7e
"""
import numpy as np
import pandas as pd
import os
import glob
import datetime as dt
import matplotlib.dates as mdate
import matplotlib.pyplot as plt
import bisect
import math
import aacgmv2
from tqdm import tqdm


def tolerant_mean(arrs):
    lens = [len(i) for i in arrs]
    arr = np.ma.empty((np.max(lens), len(arrs)))
    arr.mask = True
    for idx, l in enumerate(arrs):
        arr[:len(l), idx] = l
    return arr.mean(axis=-1), arr.std(axis=-1)



def Process_SECS(path: str, date: str, time: str, latlon: bool = True):
    '''
    This function takes in a set of Spherical Elementary Current Systems data files and processes
    them into a useful format. The processing involves extracting the datetime of each single
    file from the file name and stacking the current horizontally with the column name being 
    the longitudinal value.There are 1012 longitudinal values with only 40 unique ones. These is
    because each longitude will have values from multiple latitudes. Columns with same longitude
    are averaged over.
    Parameters
    ----------
    path : str
        Location of the SECS dataset
    date : str
        Date in str format: YYYYMMDD.
    time : str
        Starting hour in format: HH:mm.


    Returns
    -------
    None.

    '''
    
    # Validate date format
    if not date.isdigit() or len(date) != 8:
        raise ValueError("Invalid date format. Use YYYYMMDD format.")
    
    # Validate time format
    if not time[0:2].isdigit() or not time[3:5].isdigit() or len(time) != 5 or time[2] != ":":
        raise ValueError("Invalid time format. Use HH:mm format.")

    
    #Get to the location of the SECS files and define the start and stop times
    yr = date[0:4]
    month = date[4:6]
    day = date[6:8]


    path = f'{path}{yr}\{month}\{day}\*'

    time_list = []
    Current_array = np.array([])
    name = []
    
    #load the files and append into Current_array
    for file in tqdm(glob.glob(path)):
        
        N = file.find('_')#find the '_' in the filename.     
        date_time = file[N-8:N+7]

        yr = int(date_time[0:4])
        month = int(date_time[4:6])
        day = int(date_time[6:8])
        hour = int(date_time[9:11])
        minute = int(date_time[11:13])
        seconds = int(date_time[13:15])
        
        #Get the time from the file name and append to time_list
        time = dt.datetime(year=yr, month=month, day=day, hour=hour, minute=minute, second=seconds)
        time_list.append(time)
        
        #load the csv files and extract the current
        df = pd.read_csv(file, delim_whitespace=True, header=None, names=['LAT', 'LON', 'Current'])
        Cur = np.array(df.Current)
        
        
        #stack the currents together vertically
        if len(Cur) == 0:
            
            break
        else:
            if len(Current_array) == 0:
                Current_array = Cur
            else:
                Current_array = np.row_stack([Current_array, Cur])
    
    if latlon:
            
        for i in range(len(df)):
            name.append((df.LAT[i], df.LON[i]))
    else:
        for i in range(len(df)):
            name.append(df.LON[i])

    
    #Define a new DataFrame with the stacked currents, and column names being the Longitude values
    new_df = pd.DataFrame(Current_array, index=time_list, columns=name)
    
    return new_df
    



def plot_SECS(path: str, date: str, time: str, t_range: int = 6, 
              split: float = 0.3, latlon: bool = False):
    '''
    Create a keogram of SECS. This keogram is of current values averaged over latitude
    and plotted as longitude, time and current intensity. 

    Parameters
    ----------
    path : str
        path to the SECS files
    date : str
        Desired date of keogram plot in YYYYMMDD format please.
    time : str
        Time of the keogram plot in HH:mm format. 
    t_range: int 
        range of the time to be plotted (in hours). Default is 6 hours
    split: float
        How the 6 hours should be distributed around the timestep of interest.
    latlon: bool
        Whether the column name should include both the latitude and longitude or just the longitude.
        
    
    Returns
    -------
    None.

    ''' 
        
        
    current_list = []
    
    yr = date[0:4]
    month = date[4:6]
    day = date[6:8]
    
    #define the event time, which is the time of interest (ex. time of SSC or particle injection)
    event_string = f"{yr}-{month}-{day}T{time}:00.000000000"
    
    backtime = math.floor(t_range*split) #This is the difference between the start time and event time
    fronttime = t_range - backtime #This is the diff. between the end time and event time

    start = np.datetime64(event_string) - np.timedelta64(backtime, 'h') #starttime in datetime64
    end = np.datetime64(event_string) + np.timedelta64(fronttime, 'h') #end time in datetime64


    
    new_df = Process_SECS(path,date,time)
    #Group columns with the same longitude and average them all
    Grouped_mean_df = new_df.groupby(level=0, axis=1).mean()
    
    
    #Extract the x(time), y(longitude) and z(current) values for plotting the keogram
    x = Grouped_mean_df.index.values
    y = Grouped_mean_df.columns.values
    z = Grouped_mean_df.values
    
    #Normalize the current to be between 0 and 1
    normalized_z = np.where(z > 0, z, z / abs(z.min()))
    normalized_z = np.where(normalized_z <= 0, normalized_z, normalized_z / normalized_z.max())
    
    
    
    #Plot the Keogram
    start_index = bisect.bisect(x, start)
    end_index = bisect.bisect(x, end)
    time = x[start_index:end_index]
    Longitude = y
    Current = normalized_z[start_index:end_index]
    current_list.append(Current)
    Current_Array = np.array(current_list)
    Current_Avg = Current_Array.mean(axis=0)

    time, Longitude = np.meshgrid(time, Longitude)

    fig = plt.figure(figsize=(25, 15))
    ax = fig.add_subplot(111)
    img = ax.pcolormesh(time, Longitude, Current_Avg.T, vmin=-0.2, vmax=0.2, cmap='jet', shading='auto')
    cbar = plt.colorbar(img)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label('Currents [A]', fontsize=20, fontweight='bold')
    ax.set_ylabel('Longitude \n [Degrees]', fontsize=20, fontweight='bold')
    ax.set_xlabel('Time HH:MM', fontsize=20, fontweight='bold')
    ax.set_ylim([np.min(Longitude), 300])
    ax.tick_params(axis='both', which='major', labelsize=25)

    dateformat = '%H:%M'
    date_formatter = mdate.DateFormatter(dateformat)
    ax.xaxis.set_major_formatter(date_formatter)
    
    save_path = os.path.join(path,"Plots")
   
    
    if not os.path.exists(save_path):
       os.makedirs(save_path)
       print(f'Creating folder: {save_path}')

        
    filename = os.path.join(save_path,f'SECS_{date}.png')

    plt.savefig(filename)
    plt.close()

# Example of how to use the function


def SECS_MLT(path: str, date: str, time: str, save_file: bool = True):
    '''
    Define a DataFrame with the column names being latitude and longitude, and the 
    index being the time. This is done by using the Process_SECS function.
    Once defined, create a list with the row elements being [Lat, Lon, Time, Current, MLT].
    The MLT values are added on the list after using the AACGM-v2 python wrapper to calculate 
    the MLT based on time, latitude and longitude.

    Parameters
    ----------
    path : str
        Path to the SECS files 
    date : str
        date of interest in the format: YYYYMMDD
    time : str
        Time of interest in the format: HH:mm

    Returns
    -------
    None.

    '''
    
    save_path = os.path.join(path, 'Data')
    
    #Check if the file is already in save directory
    filename = f'{date}_mlt.csv'
    
    file_path = os.path.join(save_path, filename)

    # Check if the file already exists
    if os.path.exists(file_path):
        print(f"File '{filename}' already exists. It will not be overwritten.")
    else:
    
        #Get the df and define lonlat as True. This lets the column name be in the proper format
        df = Process_SECS(path, date, time, latlon = True) 
        
        combined_list = []
        
        
        # Iterate over the DataFrame and extract the elements, index, and column name
        for index, row in tqdm(df.iterrows(), total = len(df), desc = 'Calculating MLT'):
            
            mlt_ls = []
            
            for col_name, value in row.items():
                
                dtime = index
                
                #calculate mlat, mlon, and mlt from AACGMv2
                _, _, mlt = aacgmv2.get_aacgm_coord(*col_name,500,dtime, method = 'TRACE')
                
                #append all values to list
                mlt_ls.extend([mlt])
            combined_list.append(mlt_ls)
        
        # df = pd.DataFrame(data = combined_list, columns = ['time','LAT',
        #                                                    'LON','MLT',
        #                                                    'MLAT','MLON',
        #                                                    'Current'])
        
        # if save_file: df.to_csv(f'{date}_mlt.csv')
        
        return combined_list
    

        
    
if __name__ == '__main__':
    print("Can't run file as a stand alone!!!")

