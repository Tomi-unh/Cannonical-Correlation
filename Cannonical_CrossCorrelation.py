# -*- coding: utf-8 -*-
"""
Created on Thu May 18 19:32:30 2023

@author: may7e


  This script analyzes the supermag stations and creates a matrix of connections based on the correlations 
 between two stations, or nodes in this case. This is the first step in creating a network of supermag stations.
 To create the network, two sets of correlation coefficients are required. The first one being the normalization 
 correlation coefficient. This correlation coefficient is performed on the dataset over an extended period of time
 (28 days in this case). This extended correlation is used to normlaize the average degree of connection between 
 all the stations to ~ 5%. This means that at quiet times, each station (that has data) has about 5% degree 
 in the network. Normalizing the degree to 5% for all station generates station specific threshold 
 (Network_Analysis.py) that is then applied on the other set of calculated coefficients to get the time dependent
 network. 
 
 HOW TO RUN:
     Running the script is simple. 
         1.) Import main
         2.) Define the datetime string in the format: yyyymmdd-hhmm
         3.) Set all_calc to True. This both calculations in one go.
         
    ALTERNATIVELY: Just input the date time of interest in the if __name__ == '__main__' function below
                   and run the file.
 
--------------------------------------------------------------------------------------------------------------

"""
import re
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import json
from sklearn.cross_decomposition import CCA
from scipy.signal import detrend
import os
import pickle
import bisect 
from tqdm import tqdm 
from multiprocessing import Pool

json_filename = "CCA_variables.json"


#Get the global variables. This .json file only contains the path strings as
#they are the only global vriables used in the script. 
with open(json_filename, 'r') as global_file:
    global_var = json.load(global_file)
    
    

def data_cleanup(data1: pd.DataFrame, data2: pd.DataFrame, strategy: str = 'linear', inplace: bool = False) -> pd.DataFrame:
    """
    This function takes in a pandas DataFrame and cleans it up. This is done by removing values in 
    the DataFrame not necessary for use in this analysis. 
    Simple Imputer method is used for dealing with Nan values. With the default strategy being 'mean'.
    The output is just the three vector magnetic field dataset along with the time as the index
    Parameters
        ----------
        data : Pandas DataFrame to be cleaned up.
        Strategy: Interpolation is performed on the dataset to fill nans  and the strategy defines the method used.
    
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
    

    
    return data1, data2
    

    
def natural_sort_key(s):
  '''
  Key function for sorting files naturally.
  Parameters:
      - s: File name.
  Returns:
      - Key for natural sorting.
  '''

  # Split the input string into text and numeric parts
  parts = re.split(r'(\d+)', s)

  # Convert numeric parts to integers for proper numeric sorting
  parts[1::2] = map(int, parts[1::2])

  return parts
  
  
def data_generator(data, chunk_size):
  '''
    This function helps with data generation by yielding data in chunks.
    Parameters:
      ----------
        - data: Data to be chunked and yielded.
        - chunk_size: Size of each data chunk to yield.
    Yields:
      --------
        - Chunk of data.
  '''
  
  for i in range(0, len(data), chunk_size):
    yield data[i:i + chunk_size]
    
    
def fill_with_random_noise(column):
  '''
    Fill missing values in a column with random noise based on the mean and standard deviation.
    Parameters:
      ----------
        - column: A pandas Series with missing values to be filled.
  '''
  
  # Calculate the mean and standard deviation of the non-NaN values
  mean = column[~column.isna()].mean()
  std = column[~column.isna()].std()
  
  # Generate random values with the same length as the column
  random_values = np.random.normal(loc=mean, scale=std, size=len(column))
  
  # Replace NaN values with random values
  column[column.isna()] = random_values[column.isna()]
  
  
  
  
def Windowed_Correlation(df1: pd.DataFrame ,df2: pd.DataFrame, 
                         window_size: int = global_var['window_size'],
                         step: int = global_var['steps']) -> list:
    '''
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
                    

    '''
    
    
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
    #load the args 
    (main_station, compare_station, 
     path, start_time_to_timestamp, 
     stop_time_to_timestamp, Duration, 
     window_size, step)  = args
    
    #get the expected lenght of the correlation coefficients list (how many coeff. is expected)
#    print([main_station, compare_station])
    days_in_minute = 1440*Duration
    
    N = 0 #number of expected datapoints
    for i in range(0, days_in_minute - window_size + 1, step):
      N+=1
      
      
    if main_station == compare_station:
#      print(f"{main_station} and Compare Station are the same. Correlation coefficient is set to 1")
      
      corr_const = [1]*N
      
    else:
      
      primary_data = pd.read_feather(os.path.join(path, main_station))
      secondary_data = pd.read_feather(os.path.join(path, compare_station))
  
      df1 ,df2 = data_cleanup(primary_data, secondary_data)
      
      del primary_data, secondary_data  # Delete unnecessary data
      
      
      start_index1 = bisect.bisect(df1.index, start_time_to_timestamp)
      stop_index1 = bisect.bisect(df1.index, stop_time_to_timestamp)
  
      start_index2 = bisect.bisect(df2.index, start_time_to_timestamp)
      stop_index2 = bisect.bisect(df2.index, stop_time_to_timestamp)
      
      
      primary = df1.iloc[start_index1:stop_index1]
      secondary = df2.iloc[start_index2:stop_index2]
      
      if len(primary) > 0 and len(secondary) > 0: #Check is the stations have data for the time period
      
        #check for percentage of missing values in the dfs
        missing_percentages1 = (primary.isnull().sum() / len(primary)) * 100
        missing_percentages2 = (secondary.isnull().sum() / len(secondary)) * 100
        
        #if the missing percentage exceeds a %80, set the corr_const to 0
        if missing_percentages1[0] > 80 or missing_percentages2[0] > 80: 
          corr_const = [0]*N
#          print("High missing data percentage. Correlation coefficient is set to 0.")
          
        else: # fill the missing nan values with random noises. 
          
          #fill all the nans with random values that oscilates around the mean
          for col in primary.columns:
            fill_with_random_noise(primary[col])
          
          for col in secondary.columns:
            fill_with_random_noise(secondary[col])
          

          corr_const = Windowed_Correlation(primary, secondary, step = step)
        
        del primary, secondary  # Delete processed data

      else: # set the corr_const to 0 if one of the stations doesn't have any data for the time period.
#        print("No data available for the specified time period. Correlation coefficient is set to 0.")

        corr_const = [0]*N
    
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
    



def corr_matrix_parallelizer(event_datetime: str, Duration: int, path: str = global_var['path'], 
                             window_size: int = global_var['window_size'],
                             num_processes: int = 10, save_path: str = global_var['save_path'],
                             filename: str = '_CCA.pickle', chunk: int = 10000, step: int = 5):
    '''
    Imports feather files with the saved magnetometer station datasets. The data is loaded and put through 
    the windowed_correlation function and the correlation coefficients are extracted. These are store in 
    an adjecent matrix for further analysis.
    
        Parameters
        ----------
        path : Path to the feather files containing the mangetometer dataset.
        start_time: ('yyyymmdd-HHmm') Time for the start of the timeframe being analyzed.
        Duration: Duration ofthe timeframe in days. Default value is 28 days. 
        num_processes: number of proccesses used for the parallelization. Default is 10
        save_path: path to save output
        chunk: how many datapoint should the date generator used spit out at a time. Default is 10000
        filename: name of the coeff file
        step: the overlap between two consecutive windows for the CCA. Default is 5
        window_size: size of the rolling window used. Default is 128
        
        Returns
        -------
        Correlation_Matrix: Adjacent matrix (i X j) with the correlation coefficients of the ith and jth magnetometer stations.
        

    '''
    if Duration == 28:
      file_id = 'Monthly'
    else:
      file_id = 'Event_based'
#    days_to_min = Duration * 1440
    station_list = os.listdir(path)
    station_list = sorted(station_list, key = natural_sort_key)
    
    time_to_timestamp = pd.to_datetime(event_datetime, format='%Y%m%d-%H%M')
    start_time_to_timestamp = time_to_timestamp - pd.DateOffset(days = (Duration)/2)
    stop_time_to_timestamp =  time_to_timestamp + pd.DateOffset(days = (Duration)/2)
   
    #Set the path name  based on the start date
    
    start_day = event_datetime[0:8]
    save_folder = os.path.join(save_path, start_day)
    
    if not os.path.exists(save_folder):
      print(f'Creating Folder for {start_day} files' )
      os.makedirs(save_folder)
      
      
    print('Starting Canonical Correlation Analysis...')
    
    '''
    Parallelize the function. 
    '''

    
    with Pool(processes = 10) as pool:
        args_list = []
        for i, main_station in enumerate(station_list):
            for j, compare_station in enumerate(station_list):
                args_list.append((main_station, compare_station, path,
                                  start_time_to_timestamp, stop_time_to_timestamp,
                                  Duration,window_size,step))
        
        
        a = 0
        for chunk in data_generator(args_list,chunk_size = chunk):
          results = list(tqdm(pool.imap(corr_matrix, chunk), total=len(chunk), desc='Processing Item'))
          
#          save_data(results, File_SavePath)
          
          #Define the filename and path to save to. 
          File_SavePath = os.path.join(save_folder, f'{file_id}_{a}_{filename}')
          
          
          with open(File_SavePath, 'wb') as pickle_file:
              pickle.dump(results, pickle_file)
          
          del results
          
          a+=1
#    return results


def combine_pickle(datetime: str, target_phrase: str, 
                   path: str = global_var['save_path'], file_ext: str = '.pickle',
                   remove_path: bool = False, Duration: int = 28):
  
  '''
  Combines and sorts pickle files with a specific target phrase in their names into a single pickle file
  and deletes the original pickle files.
  
  Parameters:
    ----------
  -save_name: Name of the file to be saved. Defaults to combined_pickle.pkl
  -path: The directory where the pickle files are located.
  -target_phrase: The specific phrase to look for in file names.
  -file_extension: The file extension to filter by. Defaults to '.pickle'.
  -remove_path: Argument that specifies if the files should be removed after combination
                leaving only the combined file intact. Default is False. 
  
  Returns:
    --------
  -str: The path to the combined and sorted json file.
  '''
  date = datetime[0:8] #define just the date value without the time part
  
  # Create a list to hold the paths of the pickle files
  pickle_files = []
  
  save_folder = os.path.join(path, date)
  
  # Iterate over the pickle files in the directory and store their paths
  for filename in os.listdir(save_folder):
      if filename.endswith(file_ext) and target_phrase in filename:
          full_path = os.path.join(save_folder, filename)
          pickle_files.append(full_path)
  
  # Write the combined and sorted data to a single pickle file
  combined_data = []
  
  #differentiate between monthly analysis or event based analysis
  if Duration ==28:
    combined_filename = f'Month_{date}.json'
  else:
    combined_filename = f'Event_{date}.json'
    
  #Def the file path
  
  
  combined_file_path = os.path.join(save_folder, combined_filename)
  
  # Sort the list of pickle file paths
  sorted_pickle_files = sorted(pickle_files, key = natural_sort_key)
  
  print(f'Got the {file_ext} files...')
  # Read and combine the data from the sorted pickle files
  for full_path in tqdm(sorted_pickle_files, total = len(sorted_pickle_files), desc = f'Saving {file_ext} files'):
      with open(full_path, 'rb') as file:
          data = pickle.load(file)
          
          combined_data.extend(data)
          
          del data
          
      # Delete the original pickle file if specified
      if remove_path: os.remove(full_path)      
  
  with open(combined_file_path, 'w') as combined_file:
      json.dump(combined_data, combined_file)
          

  
  return_msg = f'Sucessfully Saved {combined_filename} to path: {combined_file_path}'
  
  print(return_msg)
  
  return combined_file_path
          
  
      
def main(Datetime: str, Duration: int, steps: int, Path: str = global_var['path'], 
  save_path: str = global_var['save_path'], all_calc: bool = True):
    '''
    Store the Correlaion Matrix into a pickle file and give the files the appropriate name
    
        Parameters
        ----------
        Date : Start date for the analysis of the events. Takes the form 'yyyymmdd-hhmm'.
        File_name : Name for storing the Correlation Coefficients into a pickle file once the process has been completed.
        Path: Path to the SuperMag Data location.
        all_calc: Bool value that determines if all the calculations should be done at once (monthly and day long)
                  or if they should be done individually at the user's discretion.
        
        
        Returns
        -------
        None.
    
    '''
    if all_calc == False:
        corr_matrix_parallelizer(path = Path, event_datetime = Datetime, Duration = Duration, step = steps)
        
        if Duration == 28:
          combine_pickle(datetime = Datetime, Duration = Duration, target_phrase = 'Monthly')
          
        else:
          combine_pickle(datetime = Datetime, Duration = Duration, target_phrase = 'Event')
    
    else:
        print('Calculating "Monthly" Analysis...')
        corr_matrix_parallelizer(path = Path, event_datetime = Datetime, Duration = 28, step = 5)
        
        print('Combining "Monthly" Files.')
        combine_pickle(datetime = Datetime, Duration = Duration, target_phrase = 'Monthly', remove_path = True)
        
        print('Calculating "Event" Analysis...')
        corr_matrix_parallelizer(path = Path, event_datetime = Datetime, Duration = 1, step = 2)
        
        print('Combining "Event" Files.')
        combine_pickle(datetime = Datetime, Duration = Duration, target_phrase = 'Event', remove_path = True)


if __name__ == '__main__':
    print('This script is being run as the main program...')
    main('20120101-0915', all_calc = True)

    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    
