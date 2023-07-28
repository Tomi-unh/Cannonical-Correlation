# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:27:29 2023

@author: may7e
"""



import pandas as pd
import glob 
import pickle

# os.chdir('D:')

# Data_Path = '..\\SuperMag\\*.csv' #Path to the file location for clean up.


# Data_dic = {} # Append the cleaned data into this list.

# File_SavePath = '..\\SuperMag\\'



def SuperMag_cleanUp(filename, File_SavePath):
    """
    This file is for the magnetometer station data clean up. Takes the individual station data and proccesses them 
    into a useable pandas DataFrame format. The data is stored as a dictionary into a pickle file, with the names 
    of the stations used as the key for calling each datablock.

        Parameters
        ----------
        Path : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.

    """
    
    df = pd.read_csv(filename) #load data
    df['Date_UTC'] = pd.to_datetime(df['Date_UTC'])#convert the time to datetime.
    date = df.index[-1].strftime('%Y%m%d')#Pick out the date of the year for use in filename when saving into pickle.
    
    df.rename(columns = {'Date_UTC': 'Time_Normalized'})
    
    
    """
    Open the station csv files and clean up the dataset. The time is changed from Date_UTC to datatime
    
    """
        
    df = df.set_index(['Time_Normalized']) #Set the timestamp to the index to make ordering easier.
    df.drop(['Extent','SZA','IGRF_DECL','dbz_nez','dbn_nez','dbe_nez'], axis = 1, inplace = True)  
       
    
    file_name = 'SuperMag_' + date + '.pkl' 
    Data_dic = {}
    
    # name = file[12:15] #Get the name of the station. 
    for station_name in df['IAGA'].unique():
        new_df = df.loc[df['IAGA'] == station_name] #Set the row data for each station to a new dataframe and store in a dictionary 
        
        new_df.fillna(method = 'bfill', inplace = True)
          
        Data_dic[station_name] = new_df



    with open(File_SavePath + file_name, 'wb') as pickle_file:
        pickle.dump(Data_dic, pickle_file)


        

