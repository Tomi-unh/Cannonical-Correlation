#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 18:43:28 2023

@author: tadewuyi
"""

from network_diagram import extract_corr_matrix
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import bisect 
import os



class NetworkAnalyzer:
  
  def __init__(self, data_path='../data/SuperMag', 
               global_threshold,  # User must provide a value for the global threshold
               date_str,  # User must provide the date in 'YYYYMMDD' format
               num_stations=494,  # Default value for the number of stations
               steps=5,  # Default value for the number of steps
               SM_path='../data/SuperMag'):  # Default path to SuperMag data
    
    
    self.data_path = data_path
    self.global_threshold = global_threshold
    self.date_str = date_str
    self.num_stations = num_stations
    self.steps = steps
    self.SM_path = SM_path
    
  def haversine(self, lat1: float, lon1: float, lat2: float, lon2: float) -> float: 
    """
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
    """    

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
  
  def active_stations(self, time: datetime):
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
    
    #Initialize the numbers of stations
    N_stations = 0
    for file in files:
      SM_data_path = os.path.join(self.SM_path,file)
      
      data = pd.read_feather(SM_data_path)
      
      #Get the closest time index of interest
      index = bisect.bisect_left(data['Date_UTC'],time)
      
      #Get the value of one of the magnetic field data
      value = data['dbz_geo'][index]
      
      #check to see value isn't nan and increase N_stations by 1
      if value:
        N_stations =+1
      
    return N_stations
      
      

  
  
  def Adjacent_matrix(self, correlation_result: list):
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
    #Get the correlation matrices using the extract_corr_matrix function 
    corr_dict = extract_corr_matrix(correlation_result, steps = self.steps, datetime = self.date_str)
    
    #define the adjacent matrice dictionary
    #This dictionary contains 8039 matrices each of shape 494 by 494. 8039 timesteps
    #The number of matrices can change depending on the timestep.
    adj_matrix = {}
    
    #Go through each matrix in the correlation dictionary to get the adj_matrix for each of them.
    for key, value in enumerate(corr_dict):
      
      #Subtract the global threshold from each element.
      diff = value - self.global_threshold
      
      #Stores the adjacent matrix for each correlation matrix under the same key.
      adj_matrix[key] = self.heaviside(diff)
      
    return adj_matrix
    
   
      
  def deg_connection(self, adj_matrix: dict) -> list:
    '''
    This function calculates the average degree of connections within a network that an adjacent matrix represents.
    This is based on the specfic global threshold given. This returns a list of M  (M is the total number of stations)
    values for each station with the values representing the avegerage degeree of connection of these stations 
    for a specific global threshold. 
    
    Parameter:
      ----------
      -adj_matrix: This is a dictionary that contains N matrices of shape MxM. N represnts the timestep
      while M represents the numbers of stations. 
    
    Return:
      --------
      total_array: list of M elements for each station.
                    The elements represents the average degree of connection for each station
    '''
    
    adj_len = len(adj_matrix)
    
    
    #initialize the sum total of the summation.
    total_array = []
    
    for i in range(self.num_stations):
      # 'i' select the respective station station that is being compared  to the rest.
      total = []
      for key, value in adj_matrix:
        time = datetime.strptime(key, '%Y%m%d')
        
        
        #Get the number of active stations
        N = self.active_stations(time) - 1
        
        #sum the specific row in question. These row represents the correlation 
        #value of station i with the other stations.
        sum_matrix = (sum(value[i]))/N
        
        #Append the summation into a list
        total.append(sum_matrix)
      
      #Append the average degree for each station into a list.
      total_array.append((sum(total))/adj_len)
      
    return total_array
  
  
  
  
  def total_connection(self, adj_matrix: dict) -> list:
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
    
    
    for key, value in enumerate(adj_matrix):
      
      N = self.active_stations(key)
      
      #define the total number of possible connections
      num_possible_connection  = N**2 - N
      
      tot = (np.sum(value) - np.trace(value))/num_possible_connection
      
      tot_connection_list.append(tot)
      
      
    return tot_connection_list
  
  
  
  def seperation_matrix(self):
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
    files.sort()
    
    M = len(files)
    
    #initiate an empty M x M matrix for storing the calculated values
    dist_matrix = [[0 for _ in range(M)] for _ in range(M)]
    

    for i, primary_station in enumerate(files):
      df1 = pd.read_feather(os.path.join(self.SM_path, primary_station))
      
      #get lat lon for primary station
      lat1, lon1 = df1['GEOLAT'].median(), df1['GEOLON'].median()
      
      #get the secondary station 
      for j, secondary_station in enumerate(files):
        
        df2 = pd.read_feather(os.path.join(self.SM_path, secondary_station))
        
        #get lat lon for secondary station
        lat2, lon2 = df2['GEOLAT'].median(), df2['GEOLON'].median()
        
        dist_matrix[i][j] = self.haversine(lat1,lon1,lat2,lon2) #fill the matrix with the distance calculation
        
    return dist_matrix
        
  
  
  
  
if __name__ = '__main__':
  
  '''
  Bring all the functions above together to produce various parameters of the network:
        -distance matrix 
        -average degree of connection
        -total connection 
        -adjacent matrix trace sum. This is for making sure we're analyzing properly
                  Expected value of the sum is 0.
        
  
  '''
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  