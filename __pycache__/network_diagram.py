#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 19:26:13 2023

@author: tadewuyi
"""
import networkx as nx
import numpy as np
import os
import pandas as pd
import datetime as dt
import glob
import pickle 
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
#from mpl_toolkits.basemap import Basemap
import re
import aacgmv2

 
  
  
def natural_sort_key(string):
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


def get_latlon(path: str, pattern: str = '/*.feather') -> list:
  """
  Get the latitude and longitude of each station and return a list of these values. 
  The latitdue and longitude values are sorted alphabetically for the station names. 
  
  Parameter:
    ----------
    path: Path to data folder.
    pattern: data file extension. Default is .feather
    
    Returns:
      --------
      Latitude, Longitude: 2 seperate list of latitude and longitude for the stations. 
  """
  lat = []
  lon = []
  
  station_list = []
  for file in glob.glob(f'{path}/{pattern}'):
    station_list.append(file)
  
  
  sorted_station_list = sorted(station_list, key = natural_sort_key)
  
  for file in sorted_station_list:
    df = pd.read_feather(file)
    
    
    '''
    Get the first non NAN values in the dataframe GEOLAT and GEOLON columns. This prevents putting nan values
    in the output. The Values should all be the same so it doesn't matter the specific value stored in the return
    as long as it isn't nan
    '''
    lat_index = df['GEOLAT'].first_valid_index()
    lon_index = df['GEOLON'].first_valid_index()
    
    
    
    lat.append(df.GEOLAT[lat_index])
    lon.append(df.GEOLON[lon_index])
    
  return lat, lon



def mlt_correction(SM_path):
  '''
  The MLT calculations done below to is off by some fraction. This function fixes the correction
  for each station and is then applied to the MLT calculations down below. This is done by going through
  the station dataFrame and selecting at random 10 non nan values and their respective indices from the 
  'MLT' column. The mlt calculations are made based on the selected indices and the difference between the 
  true mlt values (from the df) and the calculated mlt values are averaged to get the approximate correction.
  This assumes the difference in both values should stay relatively constant.
  
  Parameters:
    ----------
    SM_path: Path the the supermag stations files 
      
      Return:
        --------
        NONE.
  '''
  files = [file for file in os.listdir(SM_path) if '.feather' in file] #get the files in the supermag directory
  files = sorted(files, key = natural_sort_key) #sort the files in ascending order
  
  lat, lon = get_latlon(SM_path)
  
  station_diff = []
  for ind, file in tqdm(enumerate(files), total = len(files), desc = 'Calculating MLT Corrections'):
    
    df = pd.read_feather(os.path.join(SM_path, file))
    df['Date_UTC'] = pd.to_datetime(df['Date_UTC'])
    
    non_nan_indices = np.where(df['MLT'].notna())[0]
    
    random_values = np.random.choice(non_nan_indices, size=10, replace=False)
    
    mlt_ls = []
    
    for rand in random_values:
      _,_, mlt = aacgmv2.get_aacgm_coord_arr(lat[ind], lon[ind], 500, df['Date_UTC'][rand], method = 'TRACE')
      
      mlt_station = df['MLT'][rand]
      
      mlt_diff = mlt_station - mlt #get the difference between the calculated mlt and the actual mlt
      
      
      mlt_ls.append(mlt_diff)
    
    mlt_diff_avg = np.mean(mlt_ls)
    
    station_diff.append(mlt_diff_avg)
  
  diff_df = pd.DataFrame(data = station_diff, columns = ['Station diff'])
  
  diff_df.to_csv('../mlt_correction.csv')
  
  
  
    

def extract_corr_matrix(correlation_result: list, datetime: str,
                        number_of_stations: int = 494, steps: int = 5, duration: int = 2) -> dict:
    '''
    Get the correlation matrices from the results of the CCA and store these matrices into a dictionary.
    The key for the matrices are the timestep for each individual matrix. E.g, time t, t+1, t+2...
    Where time t is the beginning of the CCA analysis. 
    
    Parameters:
        ----------
        correlation_result: This is a list of lists containing the correlation coefficients between two stations.
        The list should be a square number and each sub-list in the list should represent the correlation coefficients
        between two stations.
        datetime: Date and time of the event. This is of the format 'yyyymmdd-hhmm', where the hhmm is the start
        hour and minute of the CCA analysis.
        number_of_stations: The number of stations being used for the CCA analysis. It's set to a default of 494.
        steps: The step between the first CCA analysis and the next one. Defaults to 1 for a normal rolling window.
        
        
        Returns:
        -------
        correlation_matrix: This is a dictionary of correlation matrices. Each one of these matrices corresponds 
        to a specific time step, with the appropriate label used as the key.
    '''
    

    datetime_object = dt.datetime.strptime(datetime, '%Y%m%d-%H%M')
    

    datetime_start = datetime_object - dt.timedelta(days = duration/2)
    datetime_start = datetime_start + dt.timedelta(minutes = 128)

    
    n = number_of_stations
    
    
    #Get the length of the timesteps for the correlation matrix. 
    #This is done by taking the length longest sublist in the main list.
    length_of_result = len(max(correlation_result)) 
    
    print(f'Maximum length of the row in the correlation list data is: {length_of_result}')
    
    corr_matrix = {}
    
    '''
    Make sure all the sub-lists in  the main list is all of the same length. Fill the empty ones with zeros.
    '''
    for ls in correlation_result:
        diff = length_of_result - len(ls)
        if len(ls) < length_of_result:
            ls.extend([0]*diff)

    df = pd.DataFrame(correlation_result)
    
    for col in tqdm(df.columns, total = length_of_result, desc = 'Prepping Correlation Matrix...'):
        column_data = df[col].values
        
        matrix = column_data.reshape(n,n)
        
        time_key = datetime_start + dt.timedelta(minutes = (steps)*col)
        corr_matrix[time_key] = matrix
        
        
    return corr_matrix


   
def network(adj_matrix,lat_vals,lon_vals, datetime, 
            date, fig_savepath: str = None, 
            SM_path: str = '../../../data/supermag'):

    """
  This function creates plot of networks between stations greater than 50 degrees latitude. 
  The network is drawn on a grid of latitude versus MLT. 

  Parameters:
    ----------
    adj_matrix (np.ndarray): The adjacent matrix.
    lat_vals: This is a list of the latitude of the stations used in order that they were used.
    lon_vals: List of longitude for  the stations used in the analysis in a sorted format.
    datetime (dt.datetime): The timestamp associated with the network graph.
    date: Specific date in strings of the event in question
    path (str): The path to save the generated plot. This is set to auto save to the parent directory and under
    the Network_Graph folder. If one doesn't exists, it is created. If the defined path given as an input doesn't
    exists, one is also created. 

    Returns:
    -------
        
    None
  """
    print('Loading mlt correction values')
    
    if os.path.isfile('../mlt_correction.csv'):
      mlt_cor = pd.read_csv('../mlt_correction.csv')
    else:
      mlt_correction(SM_path)
      mlt_cor = pd.read_csv('../mlt_correction.csv')
      
      
  
    _, _, mlt_vals = aacgmv2.get_aacgm_coord_arr(lat_vals, lon_vals, 500, datetime, method = 'TRACE')

    file_id = datetime.strftime('%Y%m%d-%H:%M')
    if fig_savepath is None: fig_savepath = f'../TWINS/CCA/{date}/Network_plot/Network'
    if not os.path.isdir(fig_savepath):
      print('Creating Figure Save Path')
      os.makedirs(fig_savepath)
    
    # Create a graph from the adjacency matrix
    G = nx.from_numpy_array(adj_matrix)
    
    # Assign latitude and MLT values to the nodes in the graph
    for i, node in enumerate(G.nodes()):
        # Check if the latitude value is positive
        if lat_vals[i] >= 50:
            # Add attributes to the node for latitude and MLT
            G.nodes[node]['latitude'] = lat_vals[i]
            
            
            '''
            Make sure the mlt values stay bounded in [0,24], and assign the bounded values to the nodes
            '''
            mlt_diff = mlt_vals[i] + mlt_cor['Station diff'][i]
            
            if mlt_diff <0:
              mlt_val = 24 + mlt_diff
            elif mlt_diff > 24:
              mlt_val = mlt_diff -24
            else:
              mlt_val = mlt_diff
              
            #assgin the the node
            G.nodes[node]['MLT'] = mlt_val
    
    # Filter nodes with positive latitude values
    filtered_nodes = [node for node in G.nodes() if 'latitude' in G.nodes[node] and G.nodes[node]['latitude'] >= 50]
    
    # Create a subgraph with only nodes having positive latitude
    subgraph = G.subgraph(filtered_nodes)
    
    # Plot the circular grid for latitude versus MLT
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, projection='polar')
    
    # Plot nodes on the circular grid based on latitude vs MLT
    for node in subgraph.nodes():
        # Get latitude and MLT values
        lat = G.nodes[node]['latitude']
        mlt = G.nodes[node]['MLT']
    
        # Convert latitude to radius for circular plot
        radius = lat  # Adjust the scaling if needed
    
        # Convert MLT values to radians for circular plot
        radians = mlt * 2 * 3.14159 / 24  # Convert MLT to radians
    
        # Plot nodes in polar coordinates
        ax.plot(radians, radius, marker='o', markersize=8, label=f"Node {node}")
    
    # Plot edges between connected nodes
    for edge in subgraph.edges():
        node1 = edge[0]
        node2 = edge[1]
    
        lat1 = G.nodes[node1]['latitude']
        lat2 = G.nodes[node2]['latitude']
        mlt1 = G.nodes[node1]['MLT']
        mlt2 = G.nodes[node2]['MLT']
    
        ax.plot([mlt1 * 2 * 3.14159 / 24, mlt2 * 2 * 3.14159 / 24],
                [lat1, lat2],
                color='black',
                alpha=0.5)
    
    ax.set_xticks(np.linspace(0,2*np.pi,8, endpoint = False))
    ax.set_xticklabels([f'{i * 3}:00' for i in range(8)])
    ax.set_theta_zero_location("W") #set the orientation to match that of the temperature maps
    plt.title(f'Network Graph {file_id}')
    plt.savefig(os.path.join(fig_savepath,f'Network_Graph_{file_id}.png'))
    plt.close()



def degree_hist(adj_matrix, lat: list, lon: list, 
                datetime, date: str, fig_savepath: str = None, 
                SM_path: str = '../../../data/supermag'):
  '''
  Divide the grid into section and create a histogram of the toal degree in each bin on the grid. This is a polar grid, 
  with the degree being the total number of connections made by nodes in that bin.
  '''
  
  _, _, mlt_vals = aacgmv2.get_aacgm_coord_arr(lat, lon, 500, datetime, method = 'TRACE')
  
  file_id = datetime.strftime('%Y%m%d-%H:%M')
    
    
  if fig_savepath is None: fig_savepath = f'../TWINS/CCA/{date}/Network_plot/Degree'
  
  
  if not os.path.isdir(fig_savepath):
      print('Creating Figure Save Path')
      os.makedirs(fig_savepath)
    
  # Create a graph from the adjacency matrix
  g = nx.from_numpy_array(adj_matrix)  
    
  deg_dis = dict(g.degree()) #find degree for each station
  deg = list(deg_dis.values())
    
  latitudes = np.array(lat) # Assuming 1000 random latitudes (-90 to 90 degrees)
  occurrences = np.array(deg)  # Random occurrences for each data point
  
  '''
  Define the mlt values and bound them in [0,24]
  '''
  
  print('Loading mlt correction values')
  
  if os.path.isfile('../mlt_correction.csv'):
    mlt_df = pd.read_csv('../mlt_correction.csv')
  else:
    mlt_correction(SM_path)
    mlt_df = pd.read_csv('../mlt_correction.csv')
  
  mlt_cor = mlt_df['Station diff'].to_list()
  
  mlt_values = mlt_vals + mlt_cor
  
  for ind, val in enumerate(mlt_values):
    if val < 0:
      
      ind_val = val + 24
      mlt_values[ind] = ind_val
      
    elif val > 24:
      
      ind_val = val - 24
      mlt_values[ind] = ind_val
      
    else:
      
      ind_val = val
      mlt_values[ind] = ind_val
    
  
    
  # Filtering for latitudes greater than 50
  filtered_latitudes = latitudes[latitudes > 50]
  filtered_mlt_values = mlt_values[latitudes > 50]
  filtered_occurrences = occurrences[latitudes > 50]
  
  # Creating a polar histogram
  fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))
  
  # Convert MLT to radians for plotting on a polar plot
  mlt_radians = np.radians(15 * filtered_mlt_values)  # MLT values from 0-24 hours to radians
  
  # Plotting the polar histogram based on latitude, MLT, and occurrence
  hb = ax.hist2d(mlt_radians, filtered_latitudes, bins=[24, 5], weights=filtered_occurrences, cmap='jet')
  
  ax.set_xticks(np.linspace(0,2*np.pi,8, endpoint = False))
  ax.set_xticklabels([f'{i * 3}:00' for i in range(8)])
  
  ax.set_theta_zero_location("W") #set the orientation to match that of the temperature maps 
  
  # Adding color bar to represent the occurrence frequency
  cbar = plt.colorbar(hb[3], ax=ax)
  cbar.set_label('Total Degree')
  cbar.mappable.set_clim(vmin=0, vmax=300)
  ax.set_title(f'Total Degree {file_id}')
  plt.savefig(os.path.join(fig_savepath,f'DegreeConnection_{file_id}.png'))
  plt.close()
 
  
  
  
  
def deg_centrality(adj_matrix, lat: float,lon: float, datetime, 
                   date_str: str, fig_savepath: str = None, 
                   SM_path: str = '../../../data/supermag'):
  
  
  _, _, mlt_vals = aacgmv2.get_aacgm_coord_arr(lat, lon, 500, datetime, method = 'TRACE')
  
  '''
  Define the mlt values and bound them in [0,24]
  '''
  
  print('Loading mlt correction values')
  
  if os.path.isfile('../mlt_correction.csv'):
    mlt_df = pd.read_csv('../mlt_correction.csv')
  else:
    mlt_correction(SM_path)
    mlt_df = pd.read_csv('../mlt_correction.csv')
  
  mlt_cor = mlt_df['Station diff'].to_list()
  
  mlt_values = mlt_vals + mlt_cor
  
  for ind, val in enumerate(mlt_values):
    if val < 0:
      
      ind_val = val + 24
      mlt_values[ind] = ind_val
      
    elif val > 24:
      
      ind_val = val - 24
      mlt_values[ind] = ind_val
      
    else:
      
      ind_val = val
      mlt_values[ind] = ind_val
  
  file_id = datetime.strftime('%Y%m%d-%H:%M')
  
  
  if fig_savepath is None: fig_savepath = f'../TWINS/CCA/{date_str}/Network_plot/Deg_centrality'
  
  
  if not os.path.isdir(fig_savepath):
      print('Creating Figure Save Path')
      os.makedirs(fig_savepath)
      
  lat = np.array(lat)
  
  g = nx.from_numpy_array(adj_matrix)
  
  deg_centrality = nx.degree_centrality(g)
  
  deg_vals = list(deg_centrality.values())
  
  deg_vals = np.array(deg_vals)
  
  filtered_latitudes = lat[lat > 50]
  filtered_mlt_values = mlt_values[lat > 50]
  filtered_deg = deg_vals[lat > 50]
  
  ind_val = list(enumerate(filtered_deg ))
  sorted_pairs = sorted(ind_val, key=lambda x: x[1], reverse=True)
  top_three = sorted_pairs[:3]
  
  circle_lat = [filtered_latitudes[top_three[0][0]], 
                filtered_latitudes[top_three[1][0]], 
                filtered_latitudes[top_three[2][0]]]
  
  circle_mlt = [filtered_mlt_values[top_three[0][0]], 
                filtered_mlt_values[top_three[1][0]], 
                filtered_mlt_values[top_three[2][0]]]
  
  circle_val = [top_three[0][1], top_three[1][1], top_three[2][1]]
  
  fig = plt.figure(figsize=(8, 8))
  ax = fig.add_subplot(111, polar=True)
  
  # Plot circles on the polar grid
  for value, lat, mlt in zip(circle_val, circle_lat, circle_mlt):
      # Convert latitude to radius for polar plot
#      radius = np.radians(90 - lat)  # Convert latitude to radians, 90 - lat to invert the     axis
  
      # Convert MLT values to radians for polar plot
      theta = np.radians(mlt * 15)  # Convert MLT to radians, 15 degrees per MLT hour
  
      # Calculate circle size based on value
      area = np.pi * (value**2)  # Circle size based on the value
  
      # Plot circles on the polar plot (converted to Cartesian)
      ax.scatter(theta, lat, s=100000*area, alpha=0.5)
      ax.text(theta, lat, f'{value:.2f}', color='black', ha='center', va='center')
  
  ax.set_ylim(0,90) #set the latitude limit 
  
  ax.set_xticklabels([f'{i * 3}:00' for i in range(8)])
  
  ax.set_theta_zero_location("W") #set the orientation to match that of the temperature maps
  
  ax.set_title(f'Degree Centrality {file_id}')
  plt.savefig(os.path.join(fig_savepath,f'Degree_centrality_{file_id}.png'))
  plt.close()   

    
def main(path: str, datetime: str, lat: list = None,lon: list = None, file_path: str = None, number_of_station: int = 494):
  """
  Main function for this file. Takes in the data and it's corresponding start datetime, along with the path
  and creates graphs of the network between the stations for  the specific event. The data that contains a list of correlation matrix is 
  converted into a dictionary where the key is the start time for the respective 128 minute window. The resulting dictionary is returned 
  as an output of the function. This is then looped through and passed into the network function to create the network diagram. 
  
  Parameter:
    ----------
    data: This is a list of lists. Each sub-list is the correlation for two stations for each running 
    window. 
    datetime: String of the starting time of the event in the format of 'yyyymmdd-hhmm'.
    number_of_station: Number of station for the analysis.
    path: Path to the data file. This should include the name of the file with the extension (.pickle) as well. 
  
  """
  yr = datetime[0:4]
  
  if file_path is None:  
    file_path = f'../Network_Graphs/{yr}/'
    
  if not os.path.exists(file_path):
    os.makedirs(file_path)
   
 
  with open(path, 'rb') as pickle_file:
    data = pickle.load(pickle_file) 
    
  if lat is None or lon is None:
    lat,lon = get_latlon('../data/')
      
  corr_matrix_dic = extract_corr_matrix(data,datetime)
  
  for key, value in tqdm(corr_matrix_dic.items(), desc ='Creating network graphs'):
    network(value,key, lat, lon, path = file_path)


if __name__ == '__main__':
  print('This script is being run as the main program...')
  main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    