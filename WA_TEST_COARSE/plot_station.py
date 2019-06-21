#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:51:44 2019

@author: ivica
"""


from netCDF4 import Dataset, num2date, date2index
from datetime import datetime
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pytz import timezone, utc
from bunch import Bunch

def get_station(nc, id, start_rec, end_rec):
    data = Bunch()
    data.time = num2date(nc.variables['ocean_time'][start_rec:end_rec], nc.variables['ocean_time'].units)
    data.lon = nc.variables['lon'][id]
    data.lat = nc.variables['lat'][id]
    data.hs = nc.variables['HS'][start_rec:end_rec,id]
    data.t01 = nc.variables['TM01'][start_rec:end_rec,id]
    data.t02 = nc.variables['TM02'][start_rec:end_rec,id]
    data.t10 = nc.variables['TM10'][start_rec:end_rec,id]
    data.mean_number = nc.variables['KLM'][start_rec:end_rec,id]
    data.mean_length = nc.variables['WLM'][start_rec:end_rec,id]
    data.mean_direction = nc.variables['DM'][start_rec:end_rec,id]
    data.mean_directional_spreading = nc.variables['DSPR'][start_rec:end_rec,id]
    return data

def plot_hs(p0, ID):
    title = 'Station %d (%7.4f E, %7.4f N)' %(ID, p0.lon, p0.lat)
    local = timezone('Australia/Perth')
    # put all into local time 
    time = [d.replace(tzinfo=utc).astimezone(local) for d in p0.time]
    #time = p0.time
    # Plotting 
    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle(title, fontsize=14)
    ax.plot_date(time, p0.hs, 'b')
#    ax.set_ylim([-2.5, 2.5])
    ax.set_ylabel('Significant wave height (m)', color='b')
    ax.set_xlabel('Date (AWST)', color='k')
    ax.tick_params('y', colors='b')
    ax.grid(True, linestyle='-.')
    ax.tick_params(labelcolor='b', labelsize='medium', width=3)
    ax2 = ax.twinx()
    ax2.plot_date(time, p0.mean_direction, 'r.')
#    ax2.set_ylim([-2, 2])
    ax2.set_ylabel('Mean direction (deg)', color='r')
    ax2.tick_params('y', colors='r')
    ax2.grid(True, linestyle='--')
    ax2.tick_params(labelcolor='r', labelsize='medium', width=3)
    #plt.show()
    fig.autofmt_xdate()
    outfile = os.path.join(out_path, 'Station_%d_%s.png' %(ID,time[0].strftime('%Y%m%d%H')))
    plt.savefig(outfile, dpi = 144, orientation = 'landscape', transpearent = True)
    plt.close()
    
if __name__=='__main__':
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('sta_file',  default = 'station.nc', help = 'wwm model station output file')
    parser.add_argument('--out_path',  default = './frames', help = 'output folder for figures')
    parser.add_argument('--station',  default = 1, type = int,  help = 'station ID 0 - 3')
    parser.add_argument('--start_date',  default = 'None',  help = 'start date for plot in format yyyymmddHH')
    parser.add_argument('--end_date',  default = 'None',  help = 'end date for plot in format yyyymmddHH')   

    args = parser.parse_args()

    out_path = args.out_path
    ID   = int(args.station)
    infile = args.sta_file
    start_date = args.start_date
    end_date = args.end_date
    
    nc = Dataset(infile)
    t = nc.variables['ocean_time']
    if start_date == 'None': 
        start_rec = 0
    else:
        start_rec = date2index(datetime.strptime(start_date,'%Y%m%d%H'), t, select='nearest')
    if end_date == 'None':
        end_rec = -1
    else:
        end_rec = date2index(datetime.strptime(end_date,'%Y%m%d%H'), t, select='nearest')
    print(start_rec,' -> ', end_rec)
    p0 = get_station(nc, ID, start_rec, end_rec) # station 11 (-1 in python is FPSO)
    nc.close()
    plot_hs(p0, ID)
    