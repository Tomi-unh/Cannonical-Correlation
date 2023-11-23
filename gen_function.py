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


def goes_plot(date: str, file_savepath: str, goes_path: str):
    '''
    Function to plot GOES satellite parameters like B field and the particle data. 

    Parameters
    ----------
    date : str 
        In 'yyyymmdd-hhmm' format. This is the time and date of interest.
    file_savepath : str
        Path to save the plotted figure.
    goes_path: str
        Path to the GOES file location. 

    Returns
    -------
    None.

    '''
    for sat in sat_ls:
        filename = [
        f'D:\\Random\\{sat}_magpd_19mp1_16s_{date[:8]}_{date[:8]}.nc',
        f'D:\\Random\\{sat}_magpd_19mp2_16s_{date[:8]}_{date[:8]}.nc',
        f'D:\\Random\\{sat}_magpd_19mp3_16s_{date[:8]}_{date[:8]}.nc',
        f'D:\\Random\\{sat}_magpd_19mp4_32s_{date[:8]}_{date[:8]}.nc',
        f'D:\\Random\\{sat}_magpd_19mp5_32s_{date[:8]}_{date[:8]}.nc'
        ]
        
        pattern6 = '{sat}_magneto_512ms_20161025_20161025.nc'
        
        
        
        '''
        Ignore the g at the end of every constant. It's meaningless but used to have meaning in a previous code
        and I can't be bothered to remove all of it from the code below.
        '''
        
        
        for file in filename:
            
            Ion_data = nc.Dataset(os.path.join(path,file)) #load the data 
            
            Epoch_Ion = timeconv(np.array(Ion_data.variables['time_tag'][:])/1000) #load the time data
            
            start_index = bisect.bisect_left(Epoch_Ion,start) #bisect to get the index based on time data
            
            stop_index = bisect.bisect_left(Epoch_Ion,end)  #get stop index
            
            flux_ls = []
            for i in range(1,10):
                flux_data = [np.ma.MaskedArray.filled(Ion_data.variables[f'M_{i}{val}_UDTC_UNCOR_CR'][start_index:stop_index]) for i in range(1, 10)]
                ave_flux = np.mean(flux_data, axis=0)
            
            
            Time_ls = Epoch_Ion[start_index:stop_index] #define the time range of interest 
        
        Ion_data1g = nc.Dataset(name1)
        Ion_data2g = nc.Dataset(name2)
        Ion_data3g = nc.Dataset(name3)
        Ion_data4g = nc.Dataset(name4)
        Ion_data5g = nc.Dataset(name5)
        
        
        Mag_data_g = nc.Dataset(path + pattern6)
        
        Epoch_Ion1_g = timeconv(np.array(Ion_data1g.variables['time_tag'][:])/1000)
        Epoch_Ion2_g = timeconv(np.array(Ion_data2g.variables['time_tag'][:])/1000)
        Epoch_Ion3_g = timeconv(np.array(Ion_data3g.variables['time_tag'][:])/1000)
        Epoch_Ion4_g = timeconv(np.array(Ion_data4g.variables['time_tag'][:])/1000)
        Epoch_Ion5_g = timeconv(np.array(Ion_data5g.variables['time_tag'][:])/1000)
        
        
        Epoch_Mag = timeconv(np.array(Mag_data_g.variables['time_tag'][:])/1000)
        #    Epoch_fgm = FGM_data['Epoch'][:]
        
        start_index1g = bisect.bisect_left(Epoch_Ion1_g,start)
        start_index2g = bisect.bisect_left(Epoch_Ion2_g,start)
        start_index3g = bisect.bisect_left(Epoch_Ion3_g,start)
        start_index4g = bisect.bisect_left(Epoch_Ion4_g,start)
        start_index5g = bisect.bisect_left(Epoch_Ion5_g,start)
        start_index6g = bisect.bisect_left(Epoch_Mag,start)
        
        stop_index1g = bisect.bisect_left(Epoch_Ion1_g,end)
        stop_index2g = bisect.bisect_left(Epoch_Ion2_g,end)
        stop_index3g = bisect.bisect_left(Epoch_Ion3_g,end)
        stop_index4g = bisect.bisect_left(Epoch_Ion4_g,end)
        stop_index5g = bisect.bisect_left(Epoch_Ion5_g,end)
        stop_index6g = bisect.bisect_left(Epoch_Mag,end)
        
        Epoch1_g = Epoch_Ion1_g[start_index1g:stop_index1g]
        Epoch2_g = Epoch_Ion2_g[start_index2g:stop_index2g]
        Epoch3_g = Epoch_Ion3_g[start_index3g:stop_index3g]
        Epoch4_g = Epoch_Ion4_g[start_index4g:stop_index4g]
        Epoch5_g = Epoch_Ion5_g[start_index5g:stop_index5g]
        
        Epoch6_g = Epoch_Mag[start_index6g:stop_index6g]
        
        
        # event_index6 = bisect.bisect(Epoch6_g, event)
        # event_index3 = bisect.bisect(Epoch3_g, event)
        
        
        
        P1_flux1 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_1' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux2 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_2' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux3 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_3' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux4 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_4' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux5 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_5' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux6 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_6' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux7 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_7' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux8 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_8' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        P1_flux9 = np.ma.MaskedArray.filled(Ion_data1g.variables['M_9' + val1 + '_UDTC_UNCOR_CR'][start_index1g:stop_index1g])
        
        P2_flux1 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_1' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux2 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_2' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux3 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_3' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux4 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_4' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux5 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_5' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux6 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_6' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux7 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_7' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux8 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_8' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        P2_flux9 = np.ma.MaskedArray.filled(Ion_data2g.variables['M_9' + val2 + '_UDTC_UNCOR_CR'][start_index2g:stop_index2g])
        
        P3_flux1 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_1' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux2 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_2' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux3 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_3' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux4 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_4' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux5 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_5' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux6 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_6' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux7 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_7' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux8 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_8' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        P3_flux9 = np.ma.MaskedArray.filled(Ion_data3g.variables['M_9' + val3 + '_UDTC_UNCOR_CR'][start_index3g:stop_index3g])
        
        P4_flux1 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_1' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux2 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_2' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux3 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_3' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux4 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_4' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux5 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_5' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux6 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_6' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux7 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_7' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux8 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_8' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        P4_flux9 = np.ma.MaskedArray.filled(Ion_data4g.variables['M_9' + val4 + '_UDTC_UNCOR_CR'][start_index4g:stop_index4g])
        
        P5_flux1 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_1' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux2 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_2' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux3 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_3' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux4 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_4' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux5 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_5' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux6 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_6' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux7 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_7' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux8 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_8' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        P5_flux9 = np.ma.MaskedArray.filled(Ion_data5g.variables['M_9' + val5 + '_UDTC_UNCOR_CR'][start_index5g:stop_index5g])
        
        Ave_flux1 = np.mean([P1_flux1,P1_flux2,P1_flux3,P1_flux4,P1_flux5,P1_flux6,P1_flux7,P1_flux8,P1_flux9], axis = 0)
        Ave_flux2 = np.mean([P2_flux1,P2_flux2,P2_flux3,P2_flux4,P2_flux5,P2_flux6,P2_flux7,P2_flux8,P2_flux9], axis = 0)
        Ave_flux3 = np.mean([P3_flux1,P3_flux2,P3_flux3,P3_flux4,P3_flux5,P3_flux6,P3_flux7,P3_flux8,P3_flux9], axis = 0)
        Ave_flux4 = np.mean([P4_flux1,P4_flux2,P4_flux3,P4_flux4,P4_flux5,P4_flux6,P4_flux7,P4_flux8,P4_flux9], axis = 0)
        Ave_flux5 = np.mean([P5_flux1,P5_flux2,P5_flux3,P5_flux4,P5_flux5,P5_flux6,P5_flux7,P5_flux8,P5_flux9], axis = 0)

    

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