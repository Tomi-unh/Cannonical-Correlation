# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 11:43:41 2023

@author: may7e


This is a script containig general functions for use for various things from temperature maps, to geomagnetic backgrounnd.
It's not meant to be run as a stand alone file, but imported.'
"""

from scipy.io import readsav
import numpy as np
from matplotlib import pyplot as plt
import julian
import spacepy.pycdf as cdf
import datetime as dt
import matplotlib.dates as mdate
from matplotlib.ticker import MaxNLocator
import bisect
import os
from scipy.ndimage import gaussian_filter

    

def temp_map(file_path, save_path): 
    """
    Parameters
    ----------
    file_path : STRING
        Path of the data being plotted. 
    file_save_name : STRING
        File name used for saving the generated plot

    Returns
    -------
    None.

    """
    
    """
    Define xs and ys to create the circle around the 3.3 Re region in the temperature maps. 
    """


    theta = np.linspace(0,2*np.pi, 50)

    r = 6.6

    x = r*np.cos(theta)
    y = r*np.sin(theta)


    data = readsav(file_path) # load the .sav file 
    
    temperature = data.savestruct.temparrays[0] #Stores the temperature array from the .sav file
    
    # Apply Gaussian smoothing
    smoothed_temperature = gaussian_filter(temperature, sigma=1)
    
    """
    Get the start and stop time for labelling the output file. 
    """
    time_start = julian.from_jd(data.savestruct.starttime[0], fmt = 'mjd').strftime('%H:%M:%S')
    time_stop = julian.from_jd(data.savestruct.stoptime[0], fmt = 'mjd').strftime('%H:%M:%S')
    
    date = julian.from_jd(data.savestruct.starttime[0], fmt = 'mjd').strftime('%Y%m%d')
    
    file_id = julian.from_jd(data.savestruct.stoptime[0], fmt = 'mjd').strftime('%H_%M_%S')

    fig, ax = plt.subplots(figsize=(10, 10))  
    im=ax.imshow(smoothed_temperature, cmap = 'jet', origin = 'lower', extent = (-60,20, -40, 40), vmax = 20)
    circle = plt.Circle((0,0),3, color = 'white')
    ax.add_patch(circle)
    ax.plot(x,y,'k--')
    cbar = plt.colorbar(im)
    cbar.set_label('Temperature [keV]',fontsize=15, fontweight='bold')
    cbar.ax.tick_params(labelsize=15)
    plt.xlim(-30,10)
    plt.ylim(-20,20)
    plt.xlabel('X GSM [Re]', fontsize=15, fontweight='bold')
    plt.ylabel('Y GSM [Re]', fontsize=15, fontweight='bold')
    # plt.set_ylim([np.min()])
    
    ax.tick_params(axis='both', which='major', labelsize=15, width=2, length=6)  # Adjust tick label properties
    
    # Add MLT lines
    num_mlt_sectors = 8
    mlt_interval = 2*np.pi / num_mlt_sectors
    mlt_positions = np.arange(0, 2*np.pi, mlt_interval)
    
    for mlt_position in mlt_positions:
        ax.plot([0, 40*np.cos(mlt_position)], [0, 40*np.sin(mlt_position)], 'k--', linewidth=2)
    
    plt.title(f'{time_start} - {time_stop}', fontsize= 20, fontweight='bold')
    
    file_save_name = f'ENAmap_{date}_{file_id}.png'
    
    filename = os.path.join(save_path,file_save_name)
    
    plt.savefig(filename)
    plt.close()
    

    

def plot_magnetic_indices(ax, epoch, data, ylabel, v_index, color='b', nbins=4):
    """
    Helper function to plot magnetic indices on a given subplot.

    Parameters:
    - ax: matplotlib.axes.Axes
        The subplot to plot the data on.
    - epoch: array-like
        Array of datetime objects representing the timestamps.
    - data: array-like
        Array of data values to plot.
    - ylabel: str
        Label for the y-axis.
    - v_index: int
        Index of the event timestamp in the epoch array.
    - color: str, optional
        Color of the plot (default is 'b' for blue).
    - nbins: int, optional
        Number of bins for major tick locators (default is 4).
    """
    ax.plot(epoch, data, color=color)
    ax.set_ylabel(ylabel, fontsize=20, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=20)
    
    if ylabel == 'Bz [nT]': ax.axhline(y=0, color='k', linestyle='--')  # Horizontal line at y = 0
    
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=nbins, prune='upper'))
    ax.axvline(x=epoch[v_index], color='r')  # Vertical line at the event timestamp
    ax.margins(0)



def geo_background(date: str, event_time: str, step: int = 2):
    """
    Generate a plot displaying geomagnetic indices around a specific event time.

    Parameters:
    - date: str
        Date of the event in 'yyyymmdd' format.
    - event_time: str
        Time of the event in 'hhmmss' format.
    - step: int, optional
        Time window in hours around the event time (default is 2).

    Returns:
    None
    """
    # Parse the input date and event_time
    yr, month, day = map(int, [date[:4], date[4:6], date[6:8]])
    hr, minute, seconds = map(int, [event_time[:2], event_time[2:4], event_time[4:6]])

    # Construct the path to the CDF file
    path = f'D:\\OMNI_FILES\\{date[:4]}'
    name = f'omni_hro_1min_{date[0:6]}01_v01.cdf'
    file = os.path.join(path, name)

    # Read data from the CDF file
    data = cdf.CDF(file)
    Epoch = data['Epoch'][:]
    AE = data['AE_INDEX'][:]
    Bz = data['BZ_GSM'][:]
    Pressure = data['Pressure'][:]
    Dst = data['SYM_H'][:]

    # Clean up data (replace extreme values)
    Bz[Bz >= 1000] = np.nan
    Pressure[Pressure > 30] = np.median(Pressure)

    # Calculate event time and time window
    event_time = dt.datetime(yr, month, day, hr, minute, seconds)
    start = event_time - dt.timedelta(hours=step)
    end = event_time + dt.timedelta(hours=step + 3)

    # Round start and end times to the nearest hour
    start = start.replace(minute=0, second=0, microsecond=0)
    end = end.replace(minute=0, second=0, microsecond=0)

    # Check if start month is less than event_time month
    if start.month < event_time.month:
        raise ValueError("Start time falls into the previous month. Adjust the step value.")

    # Find indices for slicing data
    start_index = bisect.bisect_left(Epoch, start)
    stop_index = bisect.bisect_left(Epoch, end)
    v_index = bisect.bisect_left(Epoch[start_index:stop_index], event_time)

    # Create the plot
    fig, (ax1, ax2, ax4) = plt.subplots(3, sharex=True, figsize=(15, 15))

    # Plot magnetic indices
    plot_magnetic_indices(ax1, Epoch[start_index:stop_index], Bz[start_index:stop_index], 'Bz [nT]', v_index)
    plot_magnetic_indices(ax2, Epoch[start_index:stop_index], AE[start_index:stop_index], 'AE Index [nT]', v_index)
    plot_magnetic_indices(ax4, Epoch[start_index:stop_index], Dst[start_index:stop_index], 'Dst Index [nT]', v_index)

    # Format x-axis with hours and minutes
    dateformat = '%H:%M'
    date_formatter = mdate.DateFormatter(dateformat)
    ax4.xaxis.set_major_formatter(date_formatter)

    # Adjust layout and display the plot
    fig.align_ylabels()
    fig.subplots_adjust(hspace=0)
    imagename = f'Geomagnetic Background {date}'
    fig.suptitle(imagename, fontsize=20, y=0.9, fontweight='bold')
    plt.show()

    # Close the CDF file
    data.close()



def plot_temperature_maps(img, path):
    """
    Parameters
    ----------
    img : list
        List of file names of the data being plotted.
    path : str
        Path where the data files are located.

    Returns
    -------
    None.

    """
    
    theta = np.linspace(0, 2*np.pi, 50)
    r = 6.6
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Create subplots with shared colorbar
    fig, axs = plt.subplots(3, 3, figsize=(15, 15), sharex=True, sharey=True)
    cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])  # Position of the colorbar

    for i, file in enumerate(img):
        data = readsav(os.path.join(path, file))  # load the .sav file 
        temperature = data.savestruct.temparrays[0]  # Stores the temperature array from the .sav file
    
        # Apply Gaussian smoothing
        smoothed_temperature = gaussian_filter(temperature, sigma=1)
    
        """
        Get the start and stop time for labelling the output file. 
        """
        time_start = julian.from_jd(data.savestruct.starttime[0], fmt='mjd').strftime('%H:%M:%S')
        time_stop = julian.from_jd(data.savestruct.stoptime[0], fmt='mjd').strftime('%H:%M:%S')
        date = julian.from_jd(data.savestruct.starttime[0], fmt='mjd').strftime('%Y%m%d')
        file_id = julian.from_jd(data.savestruct.stoptime[0], fmt='mjd').strftime('%H_%M_%S')
    
        row = i // 3
        col = i % 3
        
        axs[row, col].imshow(smoothed_temperature, cmap='jet', origin='lower', extent=(-60, 20, -40, 40), vmax=20)
        axs[row, col].add_patch(plt.Circle((0, 0), 3, color='white'))
        axs[row, col].plot(x, y, 'k--')
        axs[row, col].set_xlim(-30, 10)
        axs[row, col].set_ylim(-20, 20)
        axs[row, col].tick_params(axis='both', which='major', labelsize=10, width=2, length=6)

        # Add MLT lines
        num_mlt_sectors = 8
        mlt_interval = 2 * np.pi / num_mlt_sectors
        mlt_positions = np.arange(0, 2 * np.pi, mlt_interval)
        for mlt_position in mlt_positions:
            axs[row, col].plot([0, 40 * np.cos(mlt_position)], [0, 40 * np.sin(mlt_position)], 'k--', linewidth=1)

        # Title for each subplot
        axs[row, col].set_title(f'Subplot {i + 1}', fontsize=12, fontweight='bold')

    # Common x and y labels
    fig.text(0.5, 0.04, 'X GSM [Re]', ha='center', fontsize=12, fontweight='bold')
    fig.text(0.04, 0.5, 'Y GSM [Re]', va='center', rotation='vertical', fontsize=12, fontweight='bold')

    # Common title for the figure
    fig.suptitle(f'{time_start} - {time_stop}', fontweight='bold')

    # Colorbar
    cbar = fig.colorbar(axs[2, 2].imshow(smoothed_temperature, cmap='jet', origin='lower', extent=(-60, 20, -40, 40), vmax=20),
                        cax=cbar_ax)
    cbar.set_label('Temperature [keV]', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)

    plt.show()

if __name__ == '__main__':
    print ("Try Again. Don't Run as a stand alone file")