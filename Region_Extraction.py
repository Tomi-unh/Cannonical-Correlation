
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:28:20 2023

@author: tadewuyi
"""

import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.io import readsav
import julian
import os
import math
import pickle



def extract_features(image_data: np.ndarray, threshold: int = 3, min_size: int = 200, max_size: int = 2000) -> np.ndarray:
    '''
    This function takes in an np.ndarray image data and extracts out regions in the data that are a certain threshold
    greater than the background pixels. The threshold is a constant value multiplied by the standard deviation of the 
    image data. These regions are filtered based on the size (# of pixels in the region). A list of isolated regions of images 
    is returned.
    
    Parameters:
    ----------
    -image_data: Image array for the region extraction.
      
    Return:
    -------
      -isolated_regions: List of all the isolated regions that meets the requirement.
    '''
    
    # Calculate the standard deviation of the image
    std_dev = np.std(image_data)
    
    # Create a binary mask where values above the threshold are set to 1
    mask = np.where(image_data > (threshold)*std_dev, 1, 0).astype(np.uint8)
    
    # Use connected components labeling to identify isolated regions
    _, labels = cv2.connectedComponents(mask)
    
    isolated_regions = []
    
    for label in range(1, labels.max() + 1):
        region = (labels == label)
        region_size = np.sum(region)
    
        if min_size <= region_size <= max_size:
            # Extract the region from the original image_data
            isolated_region = image_data.copy()
            isolated_region[~region] = 0
            isolated_regions.append(isolated_region)
    
    return isolated_regions



def pieslice(img_data,angle_steps: int = 8, center_x: int = 120, center_y: int = 80):
    '''
    Parameters:
    -----------
    img_data : ndarray
        Input image data.
    angle_steps : int, optional
        Number of angle steps for pie slices. Default is 8.
    center_x : int, optional
        X-coordinate of the center. Default is 120.
    center_y : int, optional
        Y-coordinate of the center. Default is 80.

    Returns:
    --------
    selected_sections : list
        List of selected pie slice sections.
    '''
    #get the dimension of the image
    height, width = img_data.shape

    
    # Create an empty list to store the selected sections
    selected_sections = []
    
    img_mean = []
    
    #define the numbers of angles used for the pie slices
    N = int(360/angle_steps)
    
    angles = []
    
    for i in range(0,136,N):
        angles.append(i)
    
    # Create a duplicate of the angles for the other side of the image.
    angles.extend(angles)
 
    
    # Create masks for each N-degree section and apply them to the image
    for i, angle in enumerate(angles):
        # Create a new blank mask as a NumPy array
        mask = np.zeros((height, width), dtype=np.uint8)
        
        img_copy = np.copy(img_data)
        # Calculate the coordinates of the sector's bounding box
        start_angle = math.radians(angle)
        end_angle = math.radians(angle + 45)
        
        
        # Calculate the coordinates of the sector arc
        for y in range(height):
            for x in range(width):
                # Calculate the polar coordinates of the pixel relative to the image center
                dx = x - center_x
                dy = center_y - y  # Flip the y-axis direction
                
                
                if i > ((angle_steps/2) - 1):
                   dy = -dy
                
                pixel_angle = math.atan2(dy, dx)  # Calculate the angle in radians
                
                
                # Check if the pixel is within the current 45-degree section
                if start_angle <= pixel_angle < end_angle:
                    mask[y, x] = 1  # Set the pixel to white (255)
    
        # Apply the mask to the heat map to select the section
        img_copy[mask==0] = 0
        # Append the selected section to the list
        selected_sections.append(img_copy)
        
        xx = np.copy(img_copy)
        
        xx[xx==0] = np.nan
        #Get the mean of the non zero values of the image
        img_mean.append(np.nanmean(xx))
    
    return img_mean, selected_sections


def plot_regions(file_path, save_path, name):
    '''
    Plot and save isolated regions from a temperature map.

    Parameters:
    -----------
    file_path : str
        Path to the data file containing the temperature map.
    save_path : str
        Path where the images of isolated regions will be saved.
    name : str
        Name for the saved image.

    Returns:
    --------
    None.
    '''
    data = readsav(file_path)
    
    image_data = data.savestruct.temparrays[0]
    
    time_start = julian.from_jd(data.savestruct.starttime[0], fmt = 'mjd').strftime('%H:%M:%S')
    time_stop = julian.from_jd(data.savestruct.stoptime[0], fmt = 'mjd').strftime('%H:%M:%S')
    
    
    isolated_regions = extract_features(image_data)

    fig = plt.figure(figsize=(15, 5))
    
    # Plot the original image
    ax_original = fig.add_subplot(1, len(isolated_regions) + 1, 1)
    ax_original.imshow(image_data, cmap='jet', origin='lower', extent=(-60, 20, -40, 40), vmax=14)
    ax_original.set_title('Original Image')
    
    # Plot the isolated regions
    for i, region in enumerate(isolated_regions):
        ax = fig.add_subplot(1, len(isolated_regions) + 1, i + 2)
        ax.imshow(region, cmap='jet', origin='lower', extent=(-60, 20, -40, 40), vmax=14)
        ax.set_title(f'Region {i + 1}')
        
        # Save the subplot as an image file (e.g., PNG)
        plt.savefig(f'region_{i + 1}.png')
    
    filename = os.path.join(save_path, name)
    plt.title(time_start + '-' + time_stop)
    plt.savefig(filename)
    
    


def time_series(path, save_path: str = '../TWINS_Project/Region_Plot/'):
    '''
    Parameters:
    -----------
    path : str
        Path to the directory containing the data files.
    save_path : str, optional
        Path where the dictionary will be saved as a pickle file. Default is '../TWINS_Project/Region_Plot/'.

    Returns:
    --------
    None
    '''
    
    
    # get the files
    img = os.listdir(path)
    img.sort()
    
    Sector0 = []
    Sector1 = []
    Sector2 = []
    Sector3 = []
    Sector4 = []
    Sector5 = []
    Sector6 = []
    Sector7 = []
    
    for file in img:
        #Open the .sav file
        data = readsav(os.path.join(path,file))
        
        time_start = julian.from_jd(data.savestruct.starttime[0], fmt = 'mjd').strftime('%H:%M:%S')
        time_stop = julian.from_jd(data.savestruct.stoptime[0], fmt = 'mjd').strftime('%H:%M:%S')
        
        datetime = file[12:29]
        date = file[12:20]
        
        #load the temperature image
        image_data = data.savestruct.temparrays[0]
        
        img_mean, img_section = pieslice(image_data)
        
        img_section.append(img_section[0])
        
        Sector0.append(img_mean[0])
        Sector1.append(img_mean[1])
        Sector2.append(img_mean[2])
        Sector3.append(img_mean[3])
        Sector4.append(img_mean[4])
        Sector5.append(img_mean[5])
        Sector6.append(img_mean[6])
        Sector7.append(img_mean[7])
        
        fig, axs = plt.subplots(3,3, figsize=(25, 25))
        
        for i, region in enumerate(img_section):
            
            if i ==8:
                row = i //3
                col = i % 3
                # ax = fig.add_subplot(1, len(img_section) + 1, i + 2)
                axs[row, col].imshow(image_data, cmap='jet', origin='lower', extent=(-60, 20, -40, 40), vmax=14)
                axs[row, col].set_title('Original Image')
                
            else:
                row = i //3
                col = i % 3
                # ax = fig.add_subplot(1, len(img_section) + 1, i + 2)
                axs[row, col].imshow(region, cmap='jet', origin='lower', extent=(-60, 20, -40, 40), vmax=14)
                axs[row, col].set_title(f'Region {i + 1}')
                
            

        folder_path = os.path.join(save_path,f'{date}')
        
        
        if not os.path.exists(folder_path):
            # If the folder doesn't exist, create it
            os.makedirs(folder_path)
            print(f"Folder '{folder_path}' created successfully.")
        else:
            print(f"Folder '{folder_path}' already exists.")
       
        
       # Save the subplot as an image file (e.g., PNG)
        # plt.savefig(f'region_{i + 1}.png')
        filename = os.path.join(folder_path, f'{datetime}.png')
        plt.title(time_start + '-' + time_stop)
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        
    # Create a dictionary to store the results
        sect_dict = {
            '9-12': Sector0,
            '6-9': Sector1,
            '3-6': Sector2,
            '0-3': Sector3,
            '12-15': Sector4,
            '15-18': Sector5,
            '18-21': Sector6,
            '21-0': Sector7
        }
    
    
        # Save the dictionary as a pickle file
        pickle_file = os.path.join(save_path, f'{date}_sect_dict.pkl')
        with open(pickle_file, 'wb') as f:
            pickle.dump(sect_dict, f)
    
        print(f'Dictionary saved as {pickle_file}')



def plot_timeseries(date, path: str ='../TWINS_Project/Region_Plot/'):
    '''
    Plot time series data and save the plot as an image.

    Parameters:
    -----------
    date : str
        The date for which the time series data should be plotted.
    path : str, optional
        The path to the directory containing the pickle file with the data. Default is '../TWINS_Project/Region_Plot/'.

    Returns:
    --------
    None
    '''
    #Open the pickle file
    with open(os.path.join(path, f'{date}_sect_dict.pkl'), 'rb') as file:
        data = pickle.load(file)
    
    #Define the path to saving the images
    folder_path = os.path.join(path, 'Time_Series')
    
    #Create the path if it doesn't exists already
    if not os.path.exists(folder_path):
        # If the folder doesn't exist, create it
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created successfully.")
    else:
        print(f"Folder '{folder_path}' already exists.")
    
    
    fig, (ax1,ax2) = plt.subplots(2, sharex = True,  figsize =(15,15))
    
    ax1.plot(data['0-3'], color = 'b', label = '0-3')
    ax1.plot(data['3-6'], color = 'k', label = '3-6')
    ax1.plot(data['6-9'], color = 'r', label = '6-9')
    ax1.plot(data['9-12'], color = 'g', label = '9-12')
    ax1.set_ylabel('Temperature [KeV]', fontsize = 20)
    ax1.margins(x=0)
    ax1.legend()
    ax1.grid()
    
    ax2.plot(data['12-15'], color = 'b', label = '12-15')
    ax2.plot(data['15-18'], color = 'k', label = '15-18')
    ax2.plot(data['18-21'], color = 'r', label = '18-21')
    ax2.plot(data['21-0'], color = 'g', label = '21-0')
    ax2.set_ylabel('Temperature [KeV]', fontsize = 20)
    ax2.margins(x=0)
    ax2.legend()
    ax2.grid()
    

    #Save the plots
    fig.subplots_adjust(hspace= 0)#gets rid of the horizontal spacing between the subplots
    filename =  f'TimeSeries_{date}'
    fig.suptitle(filename)
    image_path = os.path.join(folder_path, filename)
    plt.savefig(image_path)
    plt.close()

    
if __name__ == '__main__':
    print ("This script shouldn't be run as is. Please import some of  the functions instead")
    