# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 20:42:41 2023

@author: may7e
"""


from gen_function import temp_map
import numpy as np
import os 
from matplotlib import pyplot as plt
import re
import julian 
from scipy.io import readsav
import datetime as dt
import pandas as pd
from tqdm.auto import tqdm
from sklearn.impute import KNNImputer
import matplotlib.dates as mdates
import matplotlib.dates as mdate
import spacepy.pycdf as cdf
import bisect

class temp_aggregate():
    def __init__(self, path, extent: tuple = (-30,10,-20,20)):
      self.path = path
      self.extent = extent
        

    def natural_sort_key(self, string):
      '''
      Key function for sorting files naturally.
      Parameters:
          - s: File name.
      Returns:
          - Key for natural sorting.
      '''
      
      # Split the input string into text and numeric parts
      parts = re.split(r'(\d+)', string)
    
      # Convert numeric parts to integers for proper numeric sorting
      parts[1::2] = map(int, parts[1::2])
    
      return parts


    def temp_max(self):
        '''
        Find the maximum temperature value in all the files being plotted and use this max value 
        to normalize all the files into a range [0,1]. The normalized temperature maps are then 
        averaged to get a single value for each map. Pixels with 0 values are excluded from the 
        averaging to give a more accurate depiction of what is going on in the temperature maps.
    
        Parameters
        ----------
        path : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        '''
        
        files = [file for file in os.listdir(self.path) if '.sav' in file]
        
        files = sorted(files, key = self.natural_sort_key)
        
        max_val = []
        
        for ind, file in tqdm(enumerate(files), total = len(files), desc = 'Getting Max value'):
            data = readsav(os.path.join(self.path,file)) 
            temperature = data.savestruct.temparrays[0] #Stores the temperature array from the .sav file 
            
            image_height, image_width = np.shape(temperature)
            # Define x and y coordinate ranges corresponding to the extent
            x_range = np.linspace(-60, 20, image_width)
            y_range = np.linspace(-40, 40, image_height)
            
            # Create boolean masks for x and y ranges
            x_mask = (x_range >= self.extent[0]) & (x_range <= self.extent[1])
            y_mask = (y_range >= self.extent[2]) & (y_range <= self.extent[3])
            
            # Apply the masks to the image
            selected_pixels = temperature[y_mask][:, x_mask]
           
            max_val.append(np.max(selected_pixels))
            
        return max(max_val)
    
    
        
    def temp_timeseries(self):
        '''
        This function takes in sets of temperature data files and converts to a timeseries.
        This is done by averaging the non zero values of the temperature maps for each timestep,
        and appending the value along with its corresponding time value to a list. The generated 
        timeseries is bounded between [0,1] from the normalization performed.

        Returns
        -------
        None.

        '''
        date = self.path[-8:]
        
        files = [file for file in os.listdir(self.path) if '.sav' in file]
        
        files = sorted(files, key = self.natural_sort_key)
        
        norm_val = self.temp_max() #normalization constant
        
        if norm_val > 200: norm_val = 200
        
        # norm_val = 1
        
        time_series = []
        time_val = []
        
        for ind, file in tqdm(enumerate(files), total = len(files), desc = 'Calculating Temperature Average'):
            if ind == 0:
                data = readsav(os.path.join(self.path,file))
                
                temperature = data.savestruct.temparrays[0] #Stores the temperature array from the .sav file
                
                image_height, image_width = np.shape(temperature)
                
                '''
                Define the region to be included in the calculations. The default is x in [-30,10], y in [-20,20]
                based in the extent(-60,20,-40,40). The (160,160) pixel image is converted to the extent range defined
                previously, and the x and y range is chosen from the extent conversion
                '''
                
                # Define x and y coordinate ranges corresponding to the extent
                x_range = np.linspace(-60, 20, image_width)
                y_range = np.linspace(-40, 40, image_height)
                
                # Create boolean masks for x and y ranges
                x_mask = (x_range >= self.extent[0]) & (x_range <= self.extent[1])
                y_mask = (y_range >= self.extent[2]) & (y_range <= self.extent[3])
                
                # Apply the masks to the image
                selected_pixels = temperature[y_mask][:, x_mask]
                
                norm_temp = selected_pixels/norm_val #normalize the temperature values 
                
                avg = norm_temp[norm_temp !=0].mean() #get the average of the non zero values
                

                #convert time strings to datetime
                dt_start = julian.from_jd(data.savestruct.starttime[0], fmt = 'mjd')
                dt_stop = julian.from_jd(data.savestruct.stoptime[0], fmt = 'mjd')
                
                #define base time for difference comparision in the else loop
                base_time = dt_stop
                
                
                #append the values to the appropriate list
                time_series.append(avg)
                time_val.append(dt_stop)
                
            else:
                data = readsav(os.path.join(self.path, file))
                
                temperature = data.savestruct.temparrays[0] #Stores the temperature array from the .sav file
                
                image_height, image_width = np.shape(temperature)
                
                # Define x and y coordinate ranges corresponding to the extent
                x_range = np.linspace(-60, 20, image_width)
                y_range = np.linspace(-40, 40, image_height)
                
                # Create boolean masks for x and y ranges
                x_mask = (x_range >= self.extent[0]) & (x_range <= self.extent[1])
                y_mask = (y_range >= self.extent[2]) & (y_range <= self.extent[3])
                
                # Apply the masks to the image
                selected_pixels = temperature[y_mask][:, x_mask]
                
                norm_temp = selected_pixels/norm_val
                avg = norm_temp[norm_temp !=0].mean()#get the average of the non zero values
                
                
                dt_start = julian.from_jd(data.savestruct.starttime[0], fmt = 'mjd')
                dt_stop = julian.from_jd(data.savestruct.stoptime[0], fmt = 'mjd')
                time_diff = (dt_stop - base_time).total_seconds()
                
                base_time = dt_stop
                
                if time_diff < 300:
                    
                    time_val.append(dt_stop)
                    time_series.append(avg)
                    
                    
                elif 300 < time_diff < 500:
                    time_elem = dt_stop + dt.timedelta(seconds = 215)
                    
                    time_val.append(time_elem)
                    time_series.append(np.nan)
                    
                elif 500< time_diff < 800:
                    time_elem = dt_stop + dt.timedelta(seconds = 215)
                    
                    time_val.append(time_elem)
                    time_series.append(np.nan)
                    
                    time_elem2 = time_elem + dt.timedelta(seconds = 215)
                    
                    time_val.append(time_elem2)
                    time_series.append(np.nan)
        
        filepath = f'../Paper/Network_plot/{date}'
        
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        
        
        filename = f'{date}_avg_temp.csv' 
        
        temp_df = pd.DataFrame({'Time': time_val, 'Temperature': time_series})
        
        temp_df_filled = temp_df.copy()
        
        # Initialize KNN imputer
        knn_imputer = KNNImputer(n_neighbors=5)  # You can adjust n_neighbors as needed
        
        # Reshape the 'Temperature' column to a 2D array
        temp_array = temp_df['Temperature'].values.reshape(-1, 1)
        temp_array_filled = knn_imputer.fit_transform(temp_array)

        
        # Update the 'Temperature' column in the filled DataFrame
        temp_df_filled['Temperature'] = temp_array_filled
        
        
        
        # Plotting columns using matplotlib
        plt.plot(temp_df['Time'],temp_df_filled['Temperature'],marker='o')
        
        # Adding labels and legend
        plt.xlabel('Time')
        plt.ylabel('Average Temperature [kev]')
        
        # Formatting the x-axis with time ticks
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=1))  # Set the interval as needed
        
        # Rotating the x-axis labels for better readability
        plt.gcf().autofmt_xdate()
        
        # Display the plot
        plt.show()
        
        
        temp_df.to_csv(os.path.join(filepath,filename))
    
    def plot_temp(self):
        
        '''
        Plot the temperature maps for the folder selected.
        '''
        date = self.path[-8:]
        
        filepath = f'../Paper/Network_plot/{date}/temp_maps'
        
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        
        files = [file for file in os.listdir(self.path) if '.sav' in file]
        
        files = sorted(files, key = self.natural_sort_key)
        
        for file in tqdm(files, total = len(files), desc = 'Plotting ENA Maps'):          
            temp_map(os.path.join(self.path,file), filepath)

    
def network_geomag_plot(date: str,
                network_path: str = '../Paper/Network_plot',
                temp_path: str = 'D:/TWINS_SAV', 
                omni_path: str = 'D:/OMNI_FILES'):
    '''
    This function takes in the temperature timeseries data, the total connection timeseries 
    of the network, db/dt, along with omni data. The output is a stacked plot of temperature, total 
    connection and geomagnetic field background. (AE,AU,AL,BZ, and Sym-H)
    
    Parameters:
        date: (yyyymmdd-hhmm) Date of the event. The plot will be based on the available temperature maps data. 
              The "hhmm" portion of the string is the event timestep. Basically the time of injection etc.
    Returns
    -------
    None.

    '''
    
    event_time = dt.datetime.strptime(date, '%Y%m%d-%H%M')
    
    #load the average temperature timeseries and define the start and stop time based on it. 
    avg_temp_data = pd.read_csv(os.path.join(network_path,f'{date[:8]}/{date[:8]}_avg_temp.csv'))
    avg_temp_data['Time'] = pd.to_datetime(avg_temp_data['Time'])
    
    start = avg_temp_data['Time'].iloc[0]
    stop = avg_temp_data['Time'].iloc[-1]
    
    #load the change in magnetic field over time data
    db_data = pd.read_csv(os.path.join(network_path,f'{date[:8]}/delta_B_{date[:8]}.csv'))
    db_data['Date_UTC'] = pd.to_datetime(db_data['Date_UTC'])
    
    #define the start,stop index for the change in magnetic field over time 
    db_start = bisect.bisect(db_data['Date_UTC'],start)
    db_stop = bisect.bisect(db_data['Date_UTC'], stop)
    db_ind = bisect.bisect(db_data['Date_UTC'], event_time)
    
    #load the total connection data. 
    tot_con_data = pd.read_csv(os.path.join(network_path,f'{date[:8]}/Network_{date[:8]}.csv')) 
    tot_con_data['Time'] = pd.to_datetime(tot_con_data['Time'])
    
    #get the start and stop index for the total degree data
    tot_start = bisect.bisect(tot_con_data['Time'],start)
    tot_stop = bisect.bisect(tot_con_data['Time'],stop)
    tot_ind = bisect.bisect(tot_con_data['Time'], event_time)
    
    
    #load the omni data and load the different variables
    omni_file = f'{date[0:4]}/omni_hro_1min_{date[:6]}01_v01.cdf'
   
    path_to_omni = os.path.join(omni_path,omni_file)
    
    omni_data = cdf.CDF(path_to_omni)
    omni_epoch = omni_data['Epoch'][:]
    AE_data = omni_data['AE_INDEX'][:]
    AL_data = omni_data['AL_INDEX'][:]
    AU_data = omni_data['AU_INDEX'][:]
    
    symh_data = omni_data['SYM_H'][:]
    Bz_data = omni_data['BZ_GSM'][:]
    Bz_data[Bz_data > 1000] = np.nan
    
    omni_start = bisect.bisect(omni_epoch, start)
    omni_stop = bisect.bisect(omni_epoch, stop)
    omni_ind = bisect.bisect(omni_epoch, event_time)
    
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(15, 15), sharex = True)
    
    ax1.plot(tot_con_data['Time'].iloc[tot_start: tot_stop], 
             tot_con_data['Total Connection'].iloc[tot_start: tot_stop], 
             color = 'b', label = 'Total Connection')
    ax1.set_ylabel('Total Connection', fontsize = 10, fontweight = 'bold')
    ax1.axvline( x = tot_con_data['Time'].iloc[tot_ind-1], color = 'k')
    ax1.margins(x=0)
    ax1.grid()

    
    ax5 = ax1.twinx()
    ax5.plot(avg_temp_data['Time'], avg_temp_data['Temperature'], 
             color = 'r', marker = 'o', label = 'Average Temperature')
    ax5.set_ylabel('Average Temperature [kev]', fontsize = 10, fontweight = 'bold')
    ax5.tick_params(axis='y', colors='red')
    ax5.yaxis.label.set_color('red')
    
    
    ax2.plot(omni_epoch[omni_start:omni_stop], Bz_data[omni_start:omni_stop],
             color = 'r')
    ax2.set_ylabel('Bz GSM [nT]',  fontsize = 10, fontweight = 'bold')
    ax2.axvline( x = omni_epoch[omni_ind], color = 'k')
    ax2.axhline(y = 0, color = 'k', linestyle = '--')
    ax2.margins(x=0)
    ax2.grid()
        
    ax6 = ax2.twinx()
    ax6.plot(omni_epoch[omni_start:omni_stop], symh_data[omni_start:omni_stop],
             color = 'b')
    ax6.set_ylabel('Sym-H [nT]',fontsize = 10, fontweight = 'bold')
    ax6.tick_params(axis='y', colors='blue')
    ax6.yaxis.label.set_color('blue')
    
    
    
    
    ax3.plot(omni_epoch[omni_start:omni_stop]
             ,AE_data[omni_start:omni_stop], 
             color = 'r', label = 'AE Index')
    ax3.plot(omni_epoch[omni_start:omni_stop]
             ,AL_data[omni_start:omni_stop], 
             color = 'b', label = 'AL Index')
    ax3.plot(omni_epoch[omni_start:omni_stop]
             ,AU_data[omni_start:omni_stop], 
             color = 'g', label = 'AU Index')
    ax3.set_ylabel('Geomagnetic Indices [nT]', fontsize = 10, fontweight = 'bold')
    ax3.axvline( x = omni_epoch[omni_ind], color = 'k')
    ax3.margins(x=0)
    ax3.legend()
    ax3.grid()
    
    
    ax4.plot(db_data['Date_UTC'].iloc[db_start : db_stop],
             db_data['mid'].iloc[db_start : db_stop], 
             color = 'b', label = 'Mid Lat.')
    ax4.plot(db_data['Date_UTC'].iloc[db_start : db_stop],
             db_data['high'].iloc[db_start : db_stop], 
             color = 'r', label = 'High Lat.')
    ax4.set_ylabel('db/dt', fontsize = 10, fontweight = 'bold')
    ax4.set_xlabel('Time UT HH:MM', fontsize = 10, fontweight = 'bold')
    ax4.axvline( x = db_data['Date_UTC'].iloc[db_ind], color = 'k')
    ax4.margins(x=0)
    ax4.legend()
    ax4.grid()
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(mdate.HourLocator(interval=1)) 
        ax.tick_params(axis='x',direction='in')
    
    dateformat = '%H:%M' #format of the time axis tick labels, with seconds being ignored
    date_formatter = mdate.DateFormatter(dateformat)
    ax4.xaxis.set_major_formatter(date_formatter)
    
    fig.align_ylabels()
    fig.subplots_adjust(hspace=0)
    fig.suptitle(f'{date[:8]} Network Connection and Geomagnetic Background')
    plt.savefig(os.path.join(network_path,f'{date[:8]}/{date[:8]}_plot.png'))
    plt.show()


if __name__ == '__main__':
    
    path = 'D:/TWINS_SAV/20150601'
    
    temp_agg = temp_aggregate(path)
    
    temp_agg.plot_temp()
    temp_agg.temp_timeseries()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        