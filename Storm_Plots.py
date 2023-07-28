# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 18:44:38 2023

@author: may7e

Plot DST and AE values for dates in the storm list beyond 2011.
"""
import pandas as pd
import numpy as np
import spacepy.pycdf as cdf
import bisect
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import os

path = '..\\stormList.csv'


event_date = pd.read_csv(path, names = ['Date'])
event_date['Date'] = pd.to_datetime(event_date['Date'])



for date in event_date['Date']:
    yr = date.strftime('%Y')
    mon = date.strftime('%m')
    date_str = date.strftime('%Y%m%d_%H%M')

    name = 'omni_hro2_1min_{yr}{month}01_v01.cdf'
    if int(yr) > 2011 and int(yr) < 2018:
        omni_data = cdf.CDF('D:\\OMNI_FILES\\' + name.format(yr = yr, month = mon))
        
        time = omni_data['Epoch'][:]
        AE = omni_data['AE_INDEX'][:]
        Sym_H = omni_data['SYM_H'][:]
        
        start = date - dt.timedelta(hours = 5)
        stop = date + dt.timedelta(hours = 5)
        
        event_index = bisect.bisect(time,date)
        event_start = bisect.bisect(time, start)
        event_stop = bisect.bisect(time, stop)
        
        fig, (ax1,ax2) = plt.subplots(2, sharex = True,  figsize =(15,15))
        
        ax1.plot(time[event_start:event_stop],AE[event_start:event_stop], color = 'b')
        ax1.set_ylabel('AE Index [nT]', fontsize = 20)
        ax1.axvline(x = time[event_index], color = 'r')
        ax1.margins(x=0)
        
        ax2.plot(time[event_start:event_stop], Sym_H[event_start:event_stop], color = 'b')
        ax2.set_xlabel('Time UT HH:MM', fontsize = 20)
        ax2.axvline(x = time[event_index], color = 'r')
        ax2.margins(x=0)
        
        dateformat = '%H:%M' #format of the time axis tick labels, with seconds being ignored
        date_formatter = mdate.DateFormatter(dateformat)
        ax2.xaxis.set_major_formatter(date_formatter)
 
        
        fig.subplots_adjust(hspace= 0)#gets rid of the horizontal spacing between the subplots
        plot_directory = '..\\Storm_plots'
        filename =  'Storm_' + date_str
        fig.suptitle(filename)
        image_path = os.path.join(plot_directory, filename)
        plt.savefig(image_path)
        plt.close()
