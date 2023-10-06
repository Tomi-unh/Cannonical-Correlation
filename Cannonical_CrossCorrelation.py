# -*- coding: utf-8 -*-
"""
Created on Thu May 18 19:32:30 2023

@author: may7e
"""

import pandas as pd
import numpy as np
from sklearn.cross_decomposition import CCA
from scipy.signal import detrend
import os
import pickle
import bisect 
from tqdm import tqdm 
from multiprocessing import Pool


def data_cleanup(data1: pd.DataFrame, data2: pd.DataFrame, strategy: str = 'mean', inplace: bool = False) -> pd.DataFrame:
    """
    This function takes in a pandas DataFrame and cleans it up. This is done by removing values in 
    the DataFrame not necessary for use in this analysis. 
    Simple Imputer method is used for dealing with Nan values. With the default strategy being 'mean'.
    The output is just the three vector magnetic field dataset along with the time as the index
    Parameters
        ----------
        data : Pandas DataFrame to be cleaned up.
        Strategy: Simple Imputer is performed on the dataset and the strategy defines the method used.
    
        Returns
        -------
        NONE: If inplace = True 
        OR (If inplace = False, which is the default)
        cleaned_data: DataFrame
            Resulting DataFrame after the cleanup has been performed in the original dataset.
    """
    
    
    'Split the DataFrame into two and detrend each'
    
    data1 = data1[['Date_UTC','dbe_geo','dbz_geo','dbn_geo']]
    data2 = data2[['Date_UTC','dbe_geo','dbz_geo','dbn_geo']]
    
    data1.set_index('Date_UTC', inplace = True)
    data2.set_index('Date_UTC', inplace = True)
    
    combined_data = pd.concat([data1, data2], axis = 1)
    
    combined_data.dropna(inplace = True)
    
    
    return combined_data
    
#    imputer = SimpleImputer(strategy  = strategy)
    
    
#    if inplace:
#        data = data[['Date_UTC','dbe_geo','dbz_geo','dbn_geo']]
#        
#        data.set_index('Date_UTC', inplace = True)
#        data1 = pd.DataFrame(imputer.fit_transform(data), columns = data.columns, index = data.index)
#        
#    else:
#        data1 = data[['Date_UTC','dbe_geo','dbz_geo','dbn_geo']]
#        data1.set_index('Date_UTC', inplace = True)
#        
#        cleaned_data = pd.DataFrame(imputer.fit_transform(data1), columns = data1.columns, index = data1.index)
#        
#        return cleaned_data
    
  
  
  
def data_generator(data, chunk_size):
  '''
  This function helps with data 
  '''
  
  for i in range(0, len(data), chunk_size):
    yield data[i:i + chunk_size]
    
    

  
def Windowed_Correlation(df1: pd.DataFrame ,df2: pd.DataFrame, window_size: int = 128, step: int = 5) -> list:
    """
    This function takes in two different DataFrames of magnetometer stations, detrends the two Datasets and performs
    a cannonical cross correlation analysis (CCA) on them. This is done in windowed segements of 128 minutes, or as otherwise specified.
    The output is an array of time dependent correlation coefficient of the two DataFrames.

        Parameters
        ----------
        df1 : Contains variables from the first magnetometer station. The index is in datetime and
              must contain magnetic field vector
        df2 : Contains variables from the first magnetometer station. The index is in datetime and
              must contain magnetic field vector
        window_size: Integer of the frequency of the rolling window performed on the DataFrames
        step: Integer of the datapoints between two windows. Defaults to 5 datapoints.
        
        Returns
        -------
        coeff_list: List of correlation coefficients after the CCA has been performed.This has a length 
                    of N - window_size. 
                    Where N: is the length of the two parsed Dataset (must be the same length).
                    

    """
    coeff_list = []
    for i in range(0,len(df1) - window_size + 1, step):
        df1_window = df1.iloc[i : i+ window_size]
        df2_window = df2.iloc[i : i+ window_size]
        
        
        detrend_df1 = detrend(df1_window[['dbe_geo','dbz_geo','dbn_geo']], axis = 0)
        detrend_df2 = detrend(df2_window[['dbe_geo','dbz_geo','dbn_geo']], axis = 0)
        
        ca = CCA(max_iter = 1500)
    
        ca.fit(detrend_df1, detrend_df2) #fit the data into a model and train.
        x_c, y_c = ca.transform(detrend_df1, detrend_df2)
        coeff = np.corrcoef(x_c[:, 0], y_c[:, 0])[0][1]
        coeff_list.append(coeff)
        
    return coeff_list
  
  
    
def corr_matrix(args):
    '''
    This function takes in two station datasets and performs the Canonical Correlation using the Windowed_Correlation 
    function and returns the correlation constant for the two pair of stations. 

    Parameters 
    ----------
    args : LIST
       Parameters must be inputted as a list, and include the following: 
           main_station (feather file). File for the main station dataset.
           compare_station (feather file). File for the secondary station dataset.
           path (string). Path to the feather files.
           start_time_to_timestamp (Timestamp). Timestamp where the analysis begins.
           days_to_min (INTEGER). Number of minutes after the start Timestamp. How long the analysis lasts for in minutes.
    Returns
    -------
    corr_const : LIST
        The correlation coefficient of the CCA performed on the two stations provided. This is the first canonical ceofficient
        used, with the other one being ignored for this analysis. 

    '''
    main_station, compare_station, path, start_time_to_timestamp, days_to_min = args
    
    primary_data = pd.read_feather(os.path.join(path, main_station))
    secondary_data = pd.read_feather(os.path.join(path, compare_station))

    combined_data = data_cleanup(primary_data, secondary_data)
    start_index = bisect.bisect(combined_data.index, start_time_to_timestamp)
    stop_index = start_index + days_to_min

    primary = combined_data.iloc[start_index:stop_index, 0:3]
    secondary = combined_data.iloc[start_index:stop_index, 3:6]
    
    del primary_data, secondary_data, combined_data  # Delete unnecessary data

    corr_const = Windowed_Correlation(primary, secondary)
    
    del primary, secondary  # Delete processed data
    
    return corr_const
  
  

def save_data(data, filename):
  '''
  This function saves data into a pickle file.
  
  Parameters:
    ----------
    -data: Data to be saved into a pickle file 
    -filename: Name of the file to save the data into.
    
    Return:
      --------
      NONE.
  '''
  with open(filename, 'wb') as pickle_file:
    pickle.dump(data, pickle_file)
    

def corr_matrix_parallelizer(path: str, start_day: str, Duration: int = 28, save_interval: int = 3600,
                             num_processes: int = 10, save_path: str = '../TWINS/CCA/',
                             filename: str = 'Month_Long_CCA.pickle', chunk: int = 100000):
    """
    Imports feather files with the saved magnetometer station datasets. The data is loaded and put through 
    the windowed_correlation function and the correlation coefficients are extracted. These are store in 
    an adjecent matrix for further analysis.
    
        Parameters
        ----------
        path : Path to the feather files containing the mangetometer dataset.
        start_time: ('yyyymmdd') Time for the start of the timeframe being analyzed.
        Duration: Duration ofthe timeframe in days. Default value is 5 days. 
        Returns
        -------
        Correlation_Matrix: Adjacent matrix (i X j) with the correlation coefficients of the ith and jth magnetometer stations.
        

    """
    
    days_to_min = Duration * 1440
    station_list = os.listdir(path)
    station_list.sort()
    
    start_time_to_timestamp = pd.to_datetime(start_day, format='%Y%m%d')

    #Set the path name 
    
    File_SavePath = os.path.join(save_path, filename)
    print('Starting Canonical Correlation Analysis...')
    
    '''
    Parallelize the function. 
    '''

    
    with Pool(processes = 10) as pool:
        args_list = []
        for i, main_station in enumerate(station_list):
            for j, compare_station in enumerate(station_list):
                args_list.append((main_station, compare_station, path, start_time_to_timestamp, days_to_min))
        
        
        for chunk in data_generator(args_list,chunk_size = chunk):
          results = list(tqdm(pool.imap(corr_matrix, chunk), total=len(chunk), desc='Processing Item'))
          
          save_data(results, File_SavePath)
          
          del results

#    return results



        
        
def main(Date: str, file_name: str,Path: str = '../data/SuperMag/', 
  save_path: str = '../TWINS/CCA/', save_interval: int = 3600):
    """
    Store the Correlaion Matrix into a pickle file and give the files the appropriate name
    
        Parameters
        ----------
        Date : Start date for the analysis of the events. Takes the form 'yyyymmdd'.
        File_name : Name for storing the Correlation Coefficients into a pickle file once the process has been completed.
        Path: Path to the SuperMag Data location.
        save_interval: Interval of periodic data dump into the pickle file. Defaults to 3600 seconds (1 hour). 
        
        Returns
        -------
        None.
    
    """
    corr_matrix_parallelizer(Path, Date)
#    Data = corr_matrix_parallelizer(Path, Date)
#    File_SavePath = os.path.join(save_path, file_name)
    
#    for data in data_generator(Data,10000):
#
#      save_data(data, File_SavePath)


if __name__ == '__main__':
    print('This script is being run as the main program...')
#    main('20120101', 'Month_Long_CCA.pkl')
    corr_matrix_parallelizer(path = '../data/SuperMag/', start_day = '20120101')
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    