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



def heaviside(array: np.ndarray):
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

def active_stations(time: datetime):
  '''
  This function gets the total number of active stations in a given time.
  
  Parameters: 
    ----------
    -time: Timestep of interest
    
    Return:
      --------
      N_stations: Total number of stations
  '''
  #Define the path to the folder
  path = '../data/SuperMag'
  
  files = os.listdir(path)
  
  #Initialize the numbers of stations
  N_stations = 0
  for file in files:
    data_path = os.path.join(path,file)
    
    data = pd.read_feather(data_path)
    
    #Get the closest time index of interest
    index = bisect.bisect_left(data['Date_UTC'],time)
    
    #Get the value of one of the magnetic field data
    value = data['dbz_geo'][index]
    
    #check to see value isn't nan and increase N_stations by 1
    if value:
      N_stations =+1
    
  return N_stations
    
    
  
def deg_connection(adj_matrix: list, global_threshold: int, date_str: str) -> int:
  '''
  This function calculates the average degree of connections within a network that an adjacent matrix represents.
  This is based on the specfic global threshold given. 
  '''
  
  
  for i in range(1, len(adj_matrix) + 1,5):
    time = datetime.strptime(date_str, '%Y%m%d') + timedelta(minutes = i)
    
    N = active_stations(time)
    
    
    

def Adjacent_matrix(correlation_result: list, date_str: str,global_threshold: int, steps: int = 5):
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
  corr_dict = extract_corr_matrix(correlation_result, steps = steps, datetime = date_str)
  
  #define the adjacent matrice dictionary
  adj_matrix = {}
  
  #Go through each matrix in the correlation dictionary to get the adj_matrix for each of them.
  for key, value in enumerate(corr_dict):
    
    #Subtract the global threshold from each element.
    diff = value - global_threshold
    
    #Stores the adjacent matrix for each correlation matrix under the same key.
    adj_matrix[key] = heaviside(diff)
    
  #Calculate 
  

    
  
  
  
  
  
  
  
  
  
  