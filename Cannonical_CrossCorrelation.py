# -*- coding: utf-8 -*-
"""
Created on Thu May 18 19:32:30 2023

@author: may7e
"""

import pandas as pd
import seaborn as sns
import numpy as np
from sklearn.cross_decomposition import CCA
from scipy.signal import detrend
import glob
import pickle


station_names =['A02', 'A04', 'A05', 'AAA', 'AAE', 'ABG', 'ABK', 'AMA', 'AMD',
                'AMS', 'AND', 'API', 'AQU', 'ARS', 'ASB', 'ASC', 'ASP', 'ATU',
                'B03', 'B07', 'B11', 'B12', 'B14', 'B18', 'B20', 'B21', 'B22',
                'B23', 'BBG', 'BCL', 'BDV', 'BEL', 'BFO', 'BJN', 'BLC', 'BMT',
                'BOR', 'BOU', 'BOX', 'BRW', 'BRZ', 'BSL', 'C01', 'C04', 'C11',
                'C12', 'CAN', 'CBB', 'CER', 'CLF', 'CMD', 'CMO', 'CNB', 'CNH',
                'CSY', 'CTA', 'CUL', 'CZT', 'DES', 'DIK', 'DMC', 'DMH', 'DOB',
                'DON', 'DOU', 'DRV', 'DRW', 'DUR', 'E05', 'EBR', 'ESK', 'EUS',
                'EYR', 'FCC', 'FHB', 'FRD', 'FRN', 'FSJ', 'FSP', 'FUR', 'FYU',
                'GAK', 'GCK', 'GDH', 'GHB', 'GHC', 'GIM', 'GNA', 'GUA', 'GUI',
                'GZH', 'HAD', 'HAN', 'HBK', 'HER', 'HLP', 'HOB', 'HON', 'HOP',
                'HRB', 'HRN', 'HUA', 'HVD', 'HYB', 'ICA', 'IGC', 'INK', 'IPM',
                'IQA', 'IRT', 'IVA', 'IZN', 'JAI', 'JYP', 'KAK', 'KAR', 'KDU',
                'KEV', 'KIL', 'KIR', 'KIV', 'KLI', 'KMH', 'KNY', 'KNZ', 'KOU',
                'KRT', 'KTB', 'KUJ', 'KUV', 'LAN', 'LCL', 'LER', 'LET', 'LIV',
                'LKW', 'LOZ', 'LRM', 'LRV', 'LVV', 'LYC', 'LYR', 'M01', 'M04',
                'M06', 'M08', 'MAB', 'MBO', 'MCQ', 'MEA', 'MEK', 'MGD', 'MIZ',
                'MMB', 'MNK', 'MOS', 'MUO', 'MUT', 'NAL', 'NAN', 'NAQ', 'NCK',
                'NEW', 'NGK', 'NOR', 'NUR', 'NVS', 'ONW', 'OSO', 'OTT', 'OUJ',
                'PAC', 'PAF', 'PAG', 'PAL', 'PBK', 'PEG', 'PEL', 'PET', 'PG1',
                'PGC', 'PHU', 'PKR', 'PPT', 'PST', 'PTK', 'PTN', 'RES', 'RIK',
                'ROE', 'SBA', 'SCO', 'SFS', 'SHE', 'SHU', 'SIT', 'SJG', 'SKT',
                'SOD', 'SOL', 'SOR', 'SPT', 'STF', 'STJ', 'SUA', 'SUW', 'SVS',
                'T03', 'T15', 'T16', 'T24', 'T25', 'T29', 'T31', 'T32', 'T33',
                'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T43', 'T58',
                'T61', 'TAL', 'TAM', 'TAR', 'TDC', 'THL', 'THY', 'TIK', 'TIR',
                'TRO', 'TRW', 'TSU', 'TUC', 'TWN', 'UMQ', 'UPN', 'UPS', 'VAL',
                'VIC', 'VRE', 'VSS', 'W02', 'WNG', 'YAK', 'YKC', 'ZAG']


def CCA_Coeff(combined_df):
    """
    This function takes in a pandas DataFrame and detrends it. The initial DataFrame is a combination of two
    seperate DataFrames that has been combined, with the first 3 columns relating 
    to the first DF and the last 3 columns to the second DF. The goal is to 
    analysis and extract the correlation between the two distinct DFs.The detrended DF is then passed into a 
    Canonical Correlation Analysis (CCA) function to return the correlation between the two Datasets.
        Parameters
        ----------
        combined_df : Pandas DataFrame
            Combination of two different SuperMag stations being compared for the CCA analysis.
    
        Returns
        -------
        ceoff: FLOAT
            First cannonical correlation coefficient of the CCA result.

    """
    
    
    'Split the DataFrame into two and detrend each'
    detrended_1 = detrend(combined_df.iloc[:,[0,1,2]], axis = 0)
    detrended_2 = detrend(combined_df.iloc[:,[3,4,5]], axis = 0)
    
    
    ca = CCA

    ca.fit(detrended_1, detrended_2) #fit the data into a model and train.
    x_c, y_c = ca.transform(detrended_1, detrended_2)
    coeff = np.corrcoef(x_c[:, 0], y_c[:, 0])[0][1]
    return coeff


def Windowed_Correlation(df1,df2, window_minute = 128):
    """
    This function takes in two different DataFrames of magnetometer stations, detrends the two Datasets and performs
    a cannonical cross correlation analysis (CCA) on them. This is done in windowed segements of 128 minutes, or as otherwise specified.
    The output is an array of time dependent correlation coefficient of the two DataFrames.

        Parameters
        ----------
        df1 : Pandas DataFrame
            Contains variables from the first magnetometer station. The index is in datetime and
            must contain magnetic field vector
        df2 : Pandas DataFrame
            Contains variables from the first magnetometer station. The index is in datetime and
            must contain magnetic field vector
        window_minute: Integer.
            Integer of the frequency of the rolling window performed on the DataFrames
        
        Returns
        -------
        Corr_coeff: Array of correlation coefficients after the CCA has been performed. 

    """
    df1.drop(['IAGA', 'GEOLON', 'GEOLAT', 'MAGON', 'MAGLAT', 'MLT', 'MCOLAT'], axis = 1, inplace = True)
    df2.drop(['IAGA', 'GEOLON', 'GEOLAT', 'MAGON', 'MAGLAT', 'MLT', 'MCOLAT'], axis = 1, inplace = True)
    
    combined = pd.concat([df1,df2], ignore_index = True, axis = 1)   
    corr_coeff = combined.rolling(window = window_minute).apply(CCA_Coeff)

    return corr_coeff
    
def corr_matrix(path):
    """
    Imports pickle files with the saved magnetometer station datasets. The data is loaded and put through 
    the windowed_correlation function and the correlation coefficients are extracted. These are store in 
    an adjecent matrix for further analysis.
    
        Parameters
        ----------
        path : String
            Path to the pickle files containing the mangetometer dataset.
        Returns
        -------
        Correlation_Matrix: Adjacent matrix (i X j) with the correlation coefficients of the ith and jth magnetometer stations.
        

    """
 
    
    # for pickle_file in glob.glob(path + '\\*.pkl'):

          
    data = pd.read_pickle(path)
    
    '''
    Get the length of time for the specific day.
    '''
    day_len = data['ABK']
    k = len(day_len)
    
    n = len(station_names)      
    
    '''Create an empty k x n x n array for creating the correlation matrix that stores the correlation coefficients
    The k dimension is the length of time in minutes. The correlation matrices useful for analysis is K - 128.
    The n dimensions are the stations used.'''
    
    corr_matrix = np.empty((k,n,n))
          
    for i in range(len(station_names)):
        primary_data = data[station_names[i]] #station data being comapred to the rest of the other stations
        
        for j in range(len(station_names)):
            secondary_data = data[station_names[j]]
            
            corr_const = Windowed_Correlation(primary_data, secondary_data) #1D array of correlation coefficients
            
            for k in range(len(corr_const)):
                corr_matrix[k][i][j] = corr_const[k]
                
    
    
    return corr_matrix 
                  
    
    
    
def main(Date, File_SavePath):
    """
    Store the Correlaion Matrix into a pickle file and give the files the appropriate name

    Parameters
    ----------
    Date : String.
        String of the date of the event being proccessed. Ex. ('01012001')
    File_SavePath : String
        Path where the file is saved

    Returns
    -------
    None.

    """
    Path = '..\\supermag\\'
    file_name = 'Corr_Matrix_' + Date
    Data = corr_matrix(Path + '\\supermag_' + Date + '.pkl')
    
    with open(File_SavePath + file_name, 'wb') as pickle_file:
        pickle.dump(Data, pickle_file)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    