# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:27:29 2023

@author: may7e
"""

"""
This file is for the magnetometer station data clean up. Takes the individual station data and proccesses them 
into a useable pandas DataFrame format. The data is stored as a dictionary into a pickle file, with the names 
of the stations used as the key for calling each datablock.
"""

import pandas as pd
import glob 
import pickle
import os

# os.chdir('D:')

Data_Path = '..\\SuperMag\\*.csv' #Path to the file location for clean up.


Data_dic = {} # Append the cleaned data into this list.

File_SavePath = 'C:/Users/may7e/Documents/GitHub/Cannonical-Correlation/'
file_name = 'SuperMag.pkl'


"""
Open the station csv files and clean up the dataset. The time is changed from Date_UTC to datatime

"""
for file in glob.glob(Data_Path):
    
    df = pd.read_csv(file) #load data
    
    name = file[12:15] #Get the name of the station. 
    
    df['Date_UTC'] = pd.to_datetime(df['Date_UTC'])
    df.drop(['Extent','IAGA','dbz_nez','dbn_nez','dbe_nez'], axis = 1, inplace = True)  

    Data_dic[name] = df
    


with open(File_SavePath + file_name, 'wb') as pickle_file:
    pickle.dump(Data_dic, pickle_file)
