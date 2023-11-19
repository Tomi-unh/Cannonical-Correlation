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
from matplotlib.patches import Circle, Arc
import julian
from matplotlib.pylab import *
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
    im=ax.imshow(smoothed_temperature, cmap = 'jet', origin = 'lower', extent = (-60,20, -40, 40), vmax = 25)
    circle = plt.Circle((0,0),3, color = 'white')
    ax.add_patch(circle)
    ax.plot(x,y,'k--')
    cbar = plt.colorbar(im)
    cbar.set_label('Temperature [keV]')
    plt.xlim(-30,10)
    plt.ylim(-20,20)
    plt.xlabel('X GSM [Re]')
    plt.ylabel('Y GSM [Re]')
    # plt.set_ylim([np.min()])
    
    # Add MLT lines
    num_mlt_sectors = 8
    mlt_interval = 2*np.pi / num_mlt_sectors
    mlt_positions = np.arange(0, 2*np.pi, mlt_interval)
    
    for mlt_position in mlt_positions:
        ax.plot([0, 40*np.cos(mlt_position)], [0, 40*np.sin(mlt_position)], 'k--')
    
    plt.title(f'{time_start} - {time_stop}')
    
    file_save_name = f'ENAmap_{date}_{file_id}.png'
    
    filename = os.path.join(save_path,file_save_name)
    
    plt.savefig(filename)
    plt.close()
    

    

def geo_background(path,plot_directory,start,end,vline):
    """
    

    Parameters
    ----------
    path : STRING
        path to the file being plotted.
    plot_directory : STRING
        Path to the folder for saving the output plot.
    start : datetime (In years, month, day, hrs, mins, secs)
        start time of the time interval of interest.
    end : datetime (In years, month, day, hrs, mins, secs)
        End time of the interval of interest
    vline : datetime (In years, month, day, hrs, mins, secs)
        Time of event where the vertical line is being drawn.

    Returns
    -------
    None.

    """
    file = path + name 
    
    data = cdf.CDF(file)
    
    Epoch = data['Epoch'][:]
    
    AE = data['AE_INDEX'][:]
    Bz = data['BZ_GSM'][:]
    Bz[Bz >= 1000] = np.median(Bz)
    Pressure = data['Pressure'][:]
    Pressure[Pressure > 30] = np.median(Pressure)
    Dst = data['SYM_H'][:]
    
    
    start_index = bisect.bisect_left(Epoch, start)
    stop_index = bisect.bisect_left(Epoch, end)
    v_index = bisect.bisect_left(Epoch, vline)
    
    
    fig, (ax1,ax2, ax4) = plt.subplots(3, sharex = True,  figsize =(15,15))
    
    
    ax1.plot(Epoch[start_index:stop_index], Bz[start_index:stop_index], color = 'b')
    ax1.set_ylabel('Bz [nT]' ,fontsize = 20, fontweight = 'bold')
    ax1.tick_params(axis = 'both',which = 'major', labelsize = 20)
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
    ax1.axhline(y = 0, color = 'k', linestyle = '--')
    ax1.axvline( x = Epoch[v_index], color = 'r')
    ax1.margins(0)
    ax1.patch.set_edgecolor('black')  
    ax1.patch.set_linewidth('5') 
    
    ax2.plot(Epoch[start_index:stop_index], AE[start_index:stop_index], color = 'b')
    ax2.set_ylabel('AE Index [nT]', fontsize = 20, fontweight = 'bold')
    ax2.tick_params(axis = 'both',which = 'major', labelsize = 20)
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
    ax2.axvline( x = Epoch[v_index], color = 'r')
    ax2.margins(0)
    ax2.patch.set_edgecolor('black')  
    ax2.patch.set_linewidth('5') 
    
    # ax3.plot(Epoch[start_index:stop_index], Pressure[start_index:stop_index], color = 'b')
    # ax3.set_ylabel('Pressure [nPa]', fontsize = 20, fontweight = 'bold')
    # ax3.tick_params(axis = 'both',which = 'major', labelsize = 20)
    # ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
    # ax3.axvline( x = Epoch[v_index], color = 'r')
    # ax3.patch.set_edgecolor('black')  
    # ax3.patch.set_linewidth('5') 
    
    ax4.plot(Epoch[start_index:stop_index], Dst[start_index:stop_index], color = 'b')
    ax4.set_ylabel('Dst Index [nT]', fontsize = 20, fontweight = 'bold')
    ax4.tick_params(axis = 'both',which = 'major', labelsize = 20)
    ax4.set_xlabel('Time', fontsize = 20, fontweight = 'bold')
    ax4.axhline(y = 0, color = 'k', linestyle = '--')
    ax4.axvline( x = Epoch[v_index], color = 'r')
    ax4.margins(0)
    ax4.patch.set_edgecolor('black')  
    ax4.patch.set_linewidth('5') 
    
    
    dateformat = '%H:%M' #format of the time axis tick labels, with seconds being ignored
    date_formatter = mdate.DateFormatter(dateformat)
    ax1.xaxis.set_major_formatter(date_formatter)
    
    fig.align_ylabels()
    fig.subplots_adjust(hspace=0)
    # plot_directory = 'C:\\Users\\may7e\\Desktop\\Event\\'
    imagename =  'Geomagnetic Background'
    filename1 =  'Geo_Background' + date
    # fig.suptitle(imagename, fontsize = 20, y = 0.9,fontweight='bold')
    image_path = os.path.join(plot_directory, filename1)
    plt.savefig(image_path)
    plt.close()




if __name__ == '__main__':
    print ("Try Again. Don't Run as a stand alone file")