#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 18:43:28 2023

@author: tadewuyi
"""

from network_diagram import (extract_corr_matrix, deg_centrality,
                             degree_hist, network, get_latlon)
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from tqdm import tqdm
import bisect 
import os
import matplotlib.pyplot as plt
import pickle
import math
import json
import re

#load the .json variable file
json_filename = 'Network_variables.json'

with open(json_filename, 'r') as global_file:
    global_var = json.load(global_file)

class NetworkAnalyzer:
  
  def __init__(self, 
               date_str: str = global_var['date_str'] ,  # User must provide the date in 'YYYYMMDD-HHMM' format
               num_stations: int = global_var['num_stations'],  # Default value for the number of stations
               num_processes: int  = global_var['num_processes'], #Default value for number of processes for multiprocessing
               steps: int = global_var['step'],  # Default value for the number of steps
               SM_path: str =global_var['SM_path']):  # Default path to SuperMag data

    self.date_str = date_str
    self.num_stations = num_stations
    self.num_processes = num_processes
    self.steps = steps
    self.SM_path = SM_path
    

    
  def cal_shortest_path(self, adj_matrix: np.ndarray, source: int):
      num_nodes = len(adj_matrix)
      dist = [math.inf] * num_nodes
      dist[source] = 0
      visited = [False] * num_nodes
  
      for _ in range(num_nodes):
          min_dist = math.inf
          min_index = -1
          for v in range(num_nodes):
              if dist[v] < min_dist and not visited[v]:
                  min_dist = dist[v]
                  min_index = v
  
          visited[min_index] = True
  
          for v in range(num_nodes):
              if not visited[v] and adj_matrix[min_index][v] > 0 and dist[v] > dist[min_index] + adj_matrix[min_index][v]:
                  dist[v] = dist[min_index] + adj_matrix[min_index][v]
  
      return dist
  
  # Calculate mean shortest path
  def mean_shortest_path(self, adj_matrix: np.ndarray) -> int:
      total_sum = 0
      num_nodes = len(adj_matrix)
  
      for i in range(num_nodes):
          distances = self.cal_shortest_path(adj_matrix, i)
          total_sum += sum(filter(lambda x: x != math.inf, distances))
  
      return total_sum / (num_nodes * (num_nodes - 1))  
    

    
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
  
  
    
  def haversine(self, lat1: float, lon1: float, lat2: float, lon2: float) -> float: 
    '''
    Calculate the great-circle distance between two points on the Earth's surface using the Haversine formula.
    
    Parameters:
    -----------
    lat1, lon1: float
        Latitude and longitude of the first point in degrees.
    lat2, lon2: float
        Latitude and longitude of the second point in degrees.
    
    Returns:
    --------
    distance: float
        The great-circle distance between the two points in kilometers.
    '''   

    # Radius of the Earth in kilometers
    R = 6371.0  # Earth's radius in kilometers

    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Differences in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Calculate the distance
    distance = R * c

    return distance
  
  
    
  def heaviside(self, array: np.ndarray):
    '''
    This is a function that performs a heaviside function on a ndarray. 
    For each element in the array, it checks whether the value is positive or negative. 
    For positive values, the function returns 1, and 0 for negative values 
    
    Parameters:
      ----------
      -array: Input array for analysis
      
      Return:
        --------
        ndarray of 0s and 1s. 
  
    '''
    return np.where(array >= 0, 1, 0)
  
  
  
  def active_stations(self, keys: list) -> list:
    '''
    This function gets the total number of active stations in a given time.
    
    Parameters: 
      ----------
      -time: Timestep of interest
      
      Return:
        --------
        N_stations: Total number of stations
    '''
    files = os.listdir(self.SM_path)
    files.sort()
    
    n_matrix = np.zeros((len(keys), len(files)))
    
    
    
    #Initialize the numbers of stations
    
    for j ,file in tqdm(enumerate(files), total = len(files), desc = 'Getting Number of Active Stations'):
      SM_data_path = os.path.join(self.SM_path,file)
      
      data = pd.read_feather(SM_data_path)
      
      #Get the closest time index of interest
     
      for i,key in enumerate(keys):
        
        index = bisect.bisect(data['Date_UTC'],key)
        #Get the value of one of the magnetic field data
        
        if index >= len(data) or index ==0:
          value = np.nan
        else:
          
          value = data['dbz_geo'][index]
        
      #check to see value isn't nan and increase N_stations by 1
        if not np.isnan(value):
          n_matrix[i][j] = 1
      
      row_ls = n_matrix.sum(axis = 1).tolist()
    return row_ls
      
      
  def delta_Bfield(self,filepath: str, start: datetime, stop: datetime):
    '''
    This function gets the changes in the magnetic field of the supermag stations as a function of time. 
    The changes in the magnetic field is over 1 minute.
    The changes are sectioned into 3 different time series
      1. The mean magnetic field magnitude change (over 1 minute) over the whole northern hemishpere stations (>50 degrees lat)
      2. The mean magnetic field magnitude change (over 1 minute) between 50-70 degrees latitude.
      3. The mean magnetic field magnitude change (over 1 minute) over 70 degrees latitude.
    '''
    
    
    if os.path.exists(filepath):
      print('Magnetic field time derivative file already exists')
    else:
      print('Calculating Magnetic field time derivative')
    
      #get the latitude values for all stations
      lat, _ = get_latlon(self.SM_path)
      
      lat = np.array(lat)
      
      #get the filenames of all stations in the supermag path
      
      files = [file for file in os.listdir(self.SM_path) if '.feather' in file]
      
      files = sorted(files, key = self.natural_sort_key)
      
      files = np.array(files) #convert files to array 
      
      '''
      fileter the latitude values 
      
      '''
      mid_condition = (lat >= 50) & (lat <= 70)
      
      #get the files that satisfy the latitude conditions
      mid_lat_station = files[mid_condition]
      high_lat_station = files[lat > 70]
      all_lat_station = files[lat >= 50]
      
      
      mid_station, high_station, all_station = [], [], []
  
      for station in tqdm(mid_lat_station, total = len(mid_lat_station), desc = 'Cal. mid lat variables'):
        df = pd.read_feather(os.path.join(self.SM_path,station))
        df['Date_UTC'] = pd.to_datetime(df['Date_UTC'])
        
        station_df = df[['Date_UTC', 'dbht']]
        
        station_df = station_df.set_index('Date_UTC')
        
        mid_station.append(station_df)
      
      print('Concatenating All Mid Lat. Station DFs')
      mid_station_concat = pd.concat(mid_station, axis = 1)
      
      print('Concatenated Mid Station DFs')
      #get the start and stop index for all data subsets
      start_ind_mid = bisect.bisect(mid_station_concat.index,start)
      stop_ind_mid = bisect.bisect(mid_station_concat.index,stop)
      
      mid_station_df = mid_station_concat.iloc[start_ind_mid:stop_ind_mid]
      
      del mid_station_concat
      
      mid_df = mid_station_df.mean(axis = 1)
      
      del mid_station_df
      
      for station in tqdm(high_lat_station, total = len(high_lat_station), desc = 'Cal. high lat variables'):
        df = pd.read_feather(os.path.join(self.SM_path,station))
        df['Date_UTC'] = pd.to_datetime(df['Date_UTC'])
        
        station_df = df[['Date_UTC', 'dbht']]
        
        station_df = station_df.set_index('Date_UTC')
        
      high_station.append(station_df)
      
      high_station_df = pd.concat(high_station, axis = 1)
      
      
      #get the start and stop index for all data subsets
      start_ind_high = bisect.bisect(high_station_df.index,start)
      stop_ind_high = bisect.bisect(high_station_df.index,stop)  
      
      high_station_df = high_station_df.iloc[start_ind_high:stop_ind_high]
      high_df = high_station_df.mean(axis = 1)
      
      del high_station_df
      
      for station in tqdm(all_lat_station, total = len(all_lat_station), desc = 'Cal. all lat variables'):
        df = pd.read_feather(os.path.join(self.SM_path,station))
        df['Date_UTC'] = pd.to_datetime(df['Date_UTC'])
        
        station_df = df[['Date_UTC', 'dbht']]
        
        station_df = station_df.set_index('Date_UTC')
        
      all_station.append(station_df)
      
      all_station_concat = pd.concat(all_station, axis = 1)
      
      
      #get the start and stop index for all data subsets
    
      start_ind_all = bisect.bisect(all_station_concat.index,start)
      stop_ind_all = bisect.bisect(all_station_concat.index,stop)
      
      all_station_df = all_station_concat.iloc[start_ind_all:stop_ind_all]
      
      del all_station_concat
      
      all_df = all_station_df.mean(axis = 1)
      
      del all_station_df
      
      #concatenate the mean values 
      concat_df = pd.concat([mid_df,high_df,all_df], axis = 1, keys = ['mid','high','all'])
      
      
      
      concat_df.to_csv(filepath)
    
    
    
    
  def Adjacent_matrix(self, corr_dict: dict, threshold: float) -> dict:
    '''
    This function takes in the CCA result which is a list of sublists. Each sublist is the CCA first cannonical 
    coefficient of the stations with respect to one specific station. This is then converted into a correlation 
    coefficient matrix for the timespan of interest. From the correlation matrix, normalized adjacent matrix is 
    obtained based on various sets of global threshold applied to each element of the corr_matrix. 
    
    Parameters:
      ----------
      -correlation_result: This is the list of sublist of the CCA analysis done on the stations.
      -date_str: This is a string of the start datetime of the CCA analysis.
      global_threshold: This is a value that sets a minimum limits what values the elements in the adjacent matrices can take.
      -steps: This is an integer of the steps between each window. In a basic rolling window, the 
              overlap between a window and the next window will be N-1, where N is the length of the 
              window. Steps defines the overlap between one window and the next N-n, where n is the step.
      
      Return:
        --------
        norm_adj_matrix: Returns the adjacent matrix as a function of time in a dictionary. 
      
    '''
    #define the adjacent matrice dictionary
    #This dictionary contains 8039 matrices each of shape 494 by 494. 8039 timesteps
    #The number of matrices can change depending on the timestep.
    adj_matrix_dict = {}
    
    #Go through each matrix in the correlation dictionary to get the adj_matrix for each of them.
    for key, value in corr_dict.items():
      
      #Subtract the global threshold from each element.
      diff = value - threshold
      
      #Stores the adjacent matrix for each correlation matrix under the same key.
      adj_matrix_dict[key] = self.heaviside(diff)
      
    return adj_matrix_dict
    
  
  
#  def parallel_deg_connection(self, adj_matrix: dict, active_stations: list) -> list:
#    '''
#    Parallelize the calculation of average degrees of connection for all stations.
#    
#    Parameters:
#      ----------
#        adj_matrix (dict): A dictionary containing N matrices of shape MxM, where N represents
#        timesteps and M represents the number of stations.
#    
#    Returns:
#      --------
#        list: A list of M elements, each representing the average degree of connection for a specific station.
#    '''
#    
#    
#    pool = multiprocessing.Pool(self.num_processes)
#    results = pool.map(self.deg_connection, [(i, adj_matrix, active_stations) for i in tqdm(range(self.num_stations))])
#    pool.close()
#    pool.join()
#    return results 
      
  def deg_connection(self, adj_matrix: dict, active_stations: list, num_station: int = 494) -> list:
    '''
    This function calculates the average degree of connections within a network that an adjacent matrix represents.
    This is based on the specfic global threshold given. This returns a list of M  (M is the total number of stations)
    values for each station with the values representing the avegerage degeree of connection of these stations 
    for a specific global threshold. 
    
    Parameter:
      ----------
      -station: This is the ith station.
      -adj_matrix: This is a dictionary that contains N matrices of shape MxM. N represnts the timestep
      while M represents the numbers of stations. 
    
    Return:
      --------
      total_array: list of M elements for each station.
                    The elements represents the average degree of connection for each station
    '''

    
    n_stations = []
    
    for station_num in range(num_station):
    
      total = []
      
      index = 0
      for key, value in adj_matrix.items():
  
        #Get the number of active stations
        N = active_stations[index] - 1
        
        index +=1
        
        if N > 0:
          #sum the specific row in question. These row represents the correlation 
          #value of station i with the other stations.
          sum_matrix = (sum(value[station_num]))/N
          
          #Append the summation into a list
          
        else:
          sum_matrix = 0 
        
        total.append(sum_matrix)
        
      tot = (sum(total))/len(adj_matrix)
      n_stations.append(tot)
      
    return n_stations
  
  
  
  def generate_datetime_list(self, starttime, num_elements, 
                             step: int = None, window_size: int = 128) ->list:
    
    '''
    Generate a list of datetime from the start time based on the numbers of element.
    The time cadence is the step.
    
    Parameter:
      ----------
      starttime: dt.datetime of the starttime of the time list generation.
      num_elements: This is the number of time indices in the time list.
      step: The step between one time element and the next. This is in minutes.
      
      Return:
        --------
        datetime_list: List of datetime values of length N, where N is the num_elements.
    '''
    if step is None: step = self.steps
    
    starttime = starttime + timedelta(minutes = window_size)
    
    datetime_list = [starttime]

    for _ in range(num_elements - 1):
        starttime += timedelta(minutes= step)
        datetime_list.append(starttime)

    return datetime_list
  
  
  
  def month_analysis(self, datetime_str: str, thresholds: tuple = None, step: int = None,
                     path: str = '../TWINS/CCA/', 
                     duration: int = 28, window_size: int = 128):
    '''
    For the CCA, a monthly analysis is required. This process calculates the CCA, and gets the 
    adjacent matrix  (based on a threshold) for the month (28 days by default). Onces generated, the adjacent matrices are 
    used for the calculation of the average degree of connection that a single station has over the 
    whole month. This average degree of connection is based on a threshold used for the calculation of 
    the adjacent matrix. This function takes in a tuple of thresholds and does the process described
    above for each threshold in tuple.
    
    Parameter:
      ----------
      datetime_str: In the format 'yyyymmdd'. This is the event date, not the start date. 
      thresholds: Tuple of the thresholds for the adjacent matrices calculations.
      step: the time step between one timestamp and the next.
      path: path to the save folder.
      duration: (Days) How long the analysis is done for. Defaults to 28 days.
      window_size: The length of the rolling window in minutes. Defaults to 128.
      
      Return:
        --------
        total_deg_connection: Total degree connection for each station as a function of the threshold used.
        
    '''
    date = datetime_str[0:8]
    
    if thresholds is None: thresholds = np.arange(0.87,1,0.01)
    
    if step is None: step = self.steps
    
    #define the start of the analysis in datetime
    start_datetime = datetime.strptime(datetime_str, '%Y%m%d-%H%M') - timedelta(days = (duration)/2)
    
    #calculate the number of datapoint in a month based on how many days used
    num_elements = 0
    duration_to_steps = duration*1440
    for i in range(0, duration_to_steps - window_size + 1, step):
      num_elements += 1
    
    datetime_list = self.generate_datetime_list(start_datetime,num_elements)
    
    
    # get the number of active stations for each timestep
    
    active_station_ls = self.active_stations(datetime_list)
    
    del datetime_list, num_elements
    #load the correlation results  for the specific date.
    
    filename = f'Month_{date}.json'
    
    file_path = os.path.join(path,date) 
    
    file = os.path.join(file_path,filename)
    
    
    #Get the number of active stations in a timestep 
    
    print('Loading Correlation Data...')
    with open(file) as f:
      data = json.load(f)
      
        #Get the correlation matrices using the extract_corr_matrix function 
    print('Calculating Correlation Matrix...')
    corr_dict = extract_corr_matrix(data, steps = self.steps, datetime = self.date_str)
    
    del data
    #Instantiate empty list for the deg. of connection for all Threshold
    total_deg_connection = []
    
    
    
    for threshold in tqdm(thresholds, total = len(thresholds), desc = 'Getting Degree of Connection...'):
      
      
      print(f'Solving for the Adjacent matrix for threshold: {threshold}')
      
      adj_dict = self.Adjacent_matrix(corr_dict, threshold = threshold)
      
      #instantiate empty list for degree connections for the specific threshold value
#      threshold_n = []
      
      print('Getting Degree of Connection...')
      
#      for station in range(self.num_stations):
      threshold_N = self.deg_connection(adj_dict, active_station_ls)
      
      
    
      total_deg_connection.append(threshold_N)
  
    tot_df = pd.DataFrame(total_deg_connection, 
                          columns =[f"Station{elem}" for elem in range(494)])
    
    
    deg_save_path = os.path.join(file_path,f'{date}_Deg_Conn.csv')
    tot_df.to_csv(deg_save_path)
    
    return total_deg_connection
    
    
    
   
    
  
  def total_connection(self, adj_matrix: dict, start_datetime: str, active_station: list,
                       duration: int = 1, window_size: int = 128) -> list:
    '''
    This takes in the adjacent matrix dictionary and returns the total number of connection 
    at each time step.
    
    Parameter: 
      ----------
      -adj_matrix: dictionary of the adjacent matrix at each timestep
      
      Return:
        --------
        tot_connection_list: list of the total connection at each timestep
    '''
    
    #initialize a list of the total connection
    tot_connection_list = []

    
    N = active_station

    index = 0
    
    
    for key, value in tqdm(adj_matrix.items(), 
                           total = len(adj_matrix), desc = 'Calculating Total Connection'):
      
      
      
      #define the total number of possible connections
      num_possible_connection  = N[index]**2 - N[index]
      
      index +=1 
      
      if num_possible_connection >0:
        tot = (np.sum(value) - np.trace(value))/num_possible_connection
        
      else:
        tot = 0
      
      
      tot_connection_list.append(tot)

    return tot_connection_list
  
  
  
  def seperation_matrix(self) -> np.ndarray:
    '''
    This function creates seperation matrix of the supermag stations.
    The elements of the matrix is the distance between sation i and j.
    
    Parameter:
      ----------
      
      Return:
        --------
        Dist_matrix: matrix of the distance between the stations
    '''
  
    files = os.listdir(self.SM_path)
    
    #sort the files to make sure the order of the matrix is correct
    files = sorted(files, key = self.natural_sort_key)
    
    M = len(files)
    
    
    location = []
    
    for file in tqdm(files, total = M, desc = 'Getting Station Location.'):
      df = pd.read_feather(os.path.join(self.SM_path,file))
      
      lat, lon = df['GEOLAT'].median(), df['GEOLON'].median()
      location.append([lat,lon])
    
    
    
    #initiate an empty M x M matrix for storing the calculated values
    dist_matrix = [[0 for _ in range(M)] for _ in range(M)]
    

    for i, latlon1 in tqdm(enumerate(location), 
                                   total = len(location), desc = 'Calculating Seperation Matrix'):
      
      
      #get lat lon for primary station
      lat1, lon1 = latlon1[0], latlon1[1]
      
      #get the secondary station 
      for j, latlon2 in enumerate(location):
        
        #get lat lon for secondary station
        lat2, lon2 = latlon2[0], latlon2[1]
        
        dist_matrix[i][j] = self.haversine(lat1,lon1,lat2,lon2) #fill the matrix with the distance calculation
    
    
    path = '../TWINS/CCA'
    name = 'distance_matrix.pickle'
    
    with open(os.path.join(path,name), 'wb') as file:
      pickle.dump(dist_matrix, file)
    
    return dist_matrix
        
  
  
  
  def avg_connection_dist(self, adj_matrix: dict, seperation_matrix: np.ndarray, active_station: list,
                          duration: int, window_size: int = 128) -> list:
    '''
    This function calculates the average connection distance between two connected stations. 
    
    Parameter:
      ----------
      adj_matrix: This is the adjacent matrix of the network.
      seperation_matrx: Matrix of the distance between two stations. 
      duration: How long the analysis is performed for in days.
      window_size: length of the rolling window in minutes. Defaults to 128 minutes.
      
      Return:
        -------
        avg_con_dist: List of the average connection distance.
    '''
    
    
    avg_con_dist = []
    
    
    N = active_station
    
    plt.plot(N)
    
    
    index = 0
    
    dis_const = np.sum(seperation_matrix) - np.trace(seperation_matrix)
    
    for key, value in tqdm(adj_matrix.items(), total = len(adj_matrix), 
                           desc = 'Calculating Average Connection Distance'):

#      n_const = N[index]**2 - N[index]
      n_const = N[index] - 1
      index +=1
      
      if n_const >0:
        
        combined = value*seperation_matrix
        combined_const = np.sum(combined) - np.trace(combined)
        
        
        avg_dist = (combined_const)/n_const
        avg_dist = (avg_dist)/dis_const

      else:
        avg_dist = 0
      
      
      avg_con_dist.append(avg_dist)
    
    return avg_con_dist
  

  
  def main(self, 
           threshold: tuple = None,
           target: int = None, 
           Duration: int = 1,
           window_size: int = 128):
    
    '''
    Bring all the functions above together to produce various parameters of the network:
        -Normalized total connection as a function of time
        -Network graph
        -Degree Histogram as a function of MLT and latitude
        -Degree Centrality for the most important stations
        
    Parameters:
      ----------
      threshold: Tuple of the thresholds to use for calculating the final adjacent matrix.
                 From the total connection, get the threshold matrix by finding the threshold coefficient
                 that gives the normalized degree of connection for each station. 
      target: The normalized degree target that each station must be close to ideally. 
      Duration: The total time lenght of the analysis. Default is 1.
      window_size: Size of the rolling window.
      
    Returns:
      --------
      Plots of the network, degree histogram, degree centrality and the normalized total degree of connection
    '''
    start_time = datetime.strptime(self.date_str, '%Y%m%d-%H%M') - timedelta(days = Duration/2)
    date = self.date_str[0:8] #used to file naming
    path = f'../TWINS/CCA/{date}' #store all the relevant files to the same location
    
    print('Getting Active Stations')
    
    
    num_elements = 0
    duration_to_steps = Duration*1440
    for i in range(0, duration_to_steps - window_size +1, self.steps):
      num_elements +=1
    
    datetime_list = self.generate_datetime_list(start_time,num_elements)
    
    N = self.active_stations(datetime_list)

    print('Calculating Average DB/DT for the Stationgs')
    
    self.delta_Bfield(filepath = os.path.join(path,f'delta_B_{date}.csv'), 
                      start = datetime_list[0], stop = datetime_list[-1])
    
    
    print('Getting the Threshold Matrix')
    
    #open and clean up the degree of connection csv file
    threshold_df = pd.read_csv(os.path.join(path,f'{date}_Deg_Conn.csv'))
    threshold_df = threshold_df.drop('Unnamed: 0', axis=1)
    
    #Define the target percent of connection. Default is 5%.
    if target is None: target = 0.05
    
    #Find in the degree of connection df, the indices of the closest deg. of Con to the target for each station.
    closest_indices_list = []
    for column_name in threshold_df.columns:
        closest_index = (threshold_df[column_name] - target).abs().idxmin()
        closest_indices_list.append(closest_index)
    
    
    #Arrray of the threshold values used to calculate the different degree of connections in the df above. 
    if threshold is None: threshold = np.arange(0.87,1,0.01) #Define one if a value isn't given.
    used_threshold = threshold
    
    threshold_matrix = np.zeros((494, 494)) #Define an empty matrix for the threshold values of each station pair.
    
    for i in range(494):
        for j in range(494):
            ind1 = closest_indices_list[i]
            ind2 = closest_indices_list[j]
            threshold1 = used_threshold[ind1]
            threshold2 = used_threshold[ind2]
            threshold_matrix[i][j] = min(threshold1, threshold2) #For the pair ij, take the lowest threshold of the two stations.
    
    '''
    Open the correlation Dictionary (or create one) and start the calculation of the Adjacent matrix for 
    every value in the correlation dictionary to form an adjacent matrix dictionary. Keep the key constant 
    between the correlation matrix and the adjacent matrix.
    '''
    #Open the correlation dictionary 
    corr_filename = f'CorrMatrix_{date}.pickle'
    corr_filepath = os.path.join(path,corr_filename)
    
    if os.path.isdir(corr_filepath):
      
      print('Loading Correlation Matrix from directory')
      with open(corr_filepath, 'rb') as file:
        
        corr_matrix_dict = pickle.load(file)
    
    else:
      
      name = f'Event_{date}.json'
      filename = os.path.join(path, name)
  
      print('Importing Correlation Data')
      with open(filename) as file:
          data = json.load(file)
  
      print('Calculating Correlation Matrix...')
      
  
      corr_matrix_dict = extract_corr_matrix(data, datetime=self.date_str, steps= self.steps, duration= Duration)
      
      del data
      
      with open(corr_filepath, 'wb') as file:
        pickle.dump(corr_matrix_dict, file)



    print('Getting Adjacent Matrix Based on Threshold Value(s)')
    adj_matrix_dict = self.Adjacent_matrix(corr_dict=corr_matrix_dict, threshold=threshold_matrix)

    time_list = list(adj_matrix_dict.keys())
#    shifted_time_list = [dt + timedelta(minutes=128) for dt in time_list]
    
    print('Calculating average db/dt')
    

    print('Getting the Total Number of Connections')
    total_connection = self.total_connection(adj_matrix_dict, start_datetime=self.date_str, active_station = N)
    tot_len = len(total_connection)
    print(f'Length of total connection data: {tot_len}')


    print('Calculating Average Connection Distance')
    
    
    dis_path = '../TWINS/CCA'
    name = 'distance_matrix.pickle'
    
    distance_matrix_path = os.path.join(dis_path,name)
    
    if os.path.isfile(distance_matrix_path):
      print('Loading Seperation Matrix')
      
      with open(distance_matrix_path, 'rb') as file:
        sep_matrix = pickle.load(file)
        
    else:
      print('Calculating Seperation Matrix.')
      
      sep_matrix = self.seperation_matrix()

    
    avg_conn_dist  = self.avg_connection_dist(adj_matrix = adj_matrix_dict, active_station = N,
                                              seperation_matrix = sep_matrix, duration = 1)

    

    print('Getting latitude and longitude values')
    
    lat, lon = get_latlon(self.SM_path)
    
    #Draw the networks

    
    for key, val in tqdm(adj_matrix_dict.items(), 
                         total = len(adj_matrix_dict), desc = 'Drawing Network Graphs...'):
      network(val, lat, lon, key, date)
    
    #Draw the degree histogram
    for key, val in tqdm(adj_matrix_dict.items(), 
                     total = len(adj_matrix_dict), desc = 'Drawing Degree Histogram...'):
      degree_hist(val, lat, lon, key, date)
    
    #Draw the degree centrality to find the most important nodes
    for key, val in tqdm(adj_matrix_dict.items(), 
                 total = len(adj_matrix_dict), desc = 'Drawing Degree Centrality...'):
      deg_centrality(val, lat, lon, key, date)



    
    network_df = pd.DataFrame({'Time': time_list, 'Connection Distance': avg_conn_dist, 
                               'Total Connection': total_connection})
  
    network_df = network_df.set_index('Time')
    
    fig, axes = plt.subplots(len(network_df.columns), 1, figsize=(8, 10), sharex=True)

    for i, column in enumerate(network_df.columns):
        axes[i].plot(network_df.index, network_df[column])
        axes[i].set_ylabel(column)

    
    plt.tight_layout()
    plt.show()
    
    network_df.to_csv(os.path.join(path, f'Network_{self.date_str[:8]}.csv'))
    
    print('DONE YAY!!!')





if __name__ == '__main__':
    analyzer = NetworkAnalyzer('19971216-2130', steps = 2)
    
    threshold = np.arange(0.87,1,0.01)
    analyzer.main(threshold) 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  