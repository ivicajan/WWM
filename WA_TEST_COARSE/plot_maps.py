#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:51:44 2019

@author: ivica
"""

from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset, num2date, date2index
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from bunch import Bunch
from pytz import timezone, utc
import matplotlib.ticker as ticker
import os

def get_history_var(fname, var, start_rec, end_rec):
    nc = Dataset(fname,'r')
    grid = Bunch()
    grid.x = nc.variables['lon'][:]
    grid.y = nc.variables['lat'][:]
    grid.e = nc.variables['ele'][:] - 1
    times = nc.variables['ocean_time']
    time = num2date(times[start_rec:end_rec], units=times.units, calendar='proleptic_gregorian')
    var = nc.variables[var][start_rec:end_rec,:]
    nc.close()
    return grid, var, time

def plot_tri_var(grid, val, trec, plotbound, clevs, cmapa, cmin, cmax, stride, title, units, orient, outfile, region):
    fig = plt.figure(figsize=(9,6), dpi=144)
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[2], urcrnrlon=plotbound[1], 
                urcrnrlat=plotbound[3], projection='merc', resolution='l', epsg='4326')
    m.drawcoastlines(linewidth=1)
    m.drawcountries()
    m.drawstates()
    m.fillcontinents(color='coral', lake_color='aqua')
    parallels = np.arange(plotbound[2], plotbound[3], np.ceil((plotbound[3]-plotbound[2])/4))
    m.drawparallels(parallels,labels = [1,0,0,0], dashes = (None, None), fontsize = 10, linewidth = 0.1)
    meridians = np.arange(plotbound[0], plotbound[1], np.ceil((plotbound[1]-plotbound[0])/4))
    m.drawmeridians(meridians, labels = [0,0,0,1], dashes = (None, None), fontsize = 10, linewidth = 0.1)
    
    cs = m.contourf(x = grid.x, y = grid.y, data = val, levels = clevs, cmap = cmapa, 
                    tri = True, latlon = True, extend='max', triangles = grid.e)
    # add colorbar.
    cbar = m.colorbar(cs, location = orient, pad = "5%", ticks = np.arange(cmin, cmax, stride), )
    cbar.set_label(units)
    # add title
    title_str = ('%s %s AWST' %(title, trec.strftime('%Y/%m/%d %H')))
    plt.title(title_str)
    plt.savefig(outfile, dpi = 144, orientation = 'landscape', transpearent = True)
    plt.close()


if __name__=='__main__':
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('history_file',  default = 'history.nc', help = 'wwm model history output file')
    parser.add_argument('--out_path',  default = './frames', help = 'output folder for figures')
    parser.add_argument('--start_date',  default = 'None',  help = 'start date for plot in format yyyymmddHH')
    parser.add_argument('--end_date',  default = 'None',  help = 'end date for plot in format yyyymmddHH')   
    parser.add_argument('--variable',  default = 'HS',  help = 'variable to plot (HS, TM01, TMO2, WLM, KLM, DM, DSPR, TPP, UBOT, ORBITAL, STOKES-BOTT-SURF-BARO-X,Y')
    parser.add_argument('--title',  default = 'Significant wave height',  help = 'title')
    parser.add_argument('--units',  default = 'm',  help = 'units for colorbar')
    parser.add_argument('--min_val',  default = 'None', help = 'min value for colorbar')
    parser.add_argument('--max_val',  default = 'None', help = 'max value for colorbar')
    parser.add_argument('--region',  default = 'wa',  help = 'region all or one of (wa, west, south, perth, geralton, shark_bay, nrc, prelude')
    parser.add_argument('--cb_location',  default = 'None',  help = 'location for colorbar (right, bottom)')
    args = parser.parse_args()

    out_path = args.out_path
    infile = args.history_file
    start_date = args.start_date
    end_date = args.end_date
    variable = args.variable
    title = args.title
    units = args.units
    region = args.region
    regions = ('wa','west','south','perth','geralton','shark_bay','nrc','prelude')
    
    #plotbound = [np.min(grid.x), np.min(grid.y), np.max(grid.x), np.max(grid.y)]
    if region=='wa':   plotbound = [108, 128, -36, -11]
    if region=='west': plotbound = [110, 116, -34.5, -21.5]
    if region=='south': plotbound = [113, 118, -36.5, -32]
    if region=='perth': plotbound = [114, 116, -33.75, -30]
    if region=='geralton': plotbound = [112, 115.5, -30, -27]
    if region=='shark_bay': plotbound = [111, 114.5, -27, -22]
    if region=='nrc': plotbound = [112.5, 120, -22.5, -17.5]
    if region=='prelude': plotbound = [121, 125, -17, -12]

    if args.cb_location == 'None': cb_location = 'right'

    start_rec = 0
    end_rec = -1
    nc = Dataset(infile)
    t = nc.variables['ocean_time']
    if start_date!='None': start_rec = date2index(datetime.strptime(start_date,'%Y%m%d%H'), t, select='nearest')
    if end_date!='None': end_rec = date2index(datetime.strptime(end_date,'%Y%m%d%H'), t, select='nearest')
    nc.close()
    print(start_rec,' -> ', end_rec)

    grid, var, time = get_history_var(infile, variable, start_rec, end_rec)
<<<<<<< HEAD
    #cmapa = make_mycolormap()
    cmapa = 'rainbow'
    if args.cb_location == 'None': cb_location = 'right'
=======
>>>>>>> refs/remotes/origin/ivica
    nrecs, npoints = np.shape(var)
    local = timezone('Australia/Perth')
    # put all into local time 
    time_local = [d.replace(tzinfo=utc).astimezone(local) for d in time]
    for i in range(0,nrecs):
        val = var[i,:]
        trec = time_local[i]
        if args.min_val == 'None': 
            vmin = np.min(val)
        else:
            vmin = np.float(args.min_val)
        if args.max_val == 'None': 
            vmax = np.max(val)
        else:
            vmax = np.float(args.max_val)    
        stride = np.ceil(np.float((vmax-vmin)/5))
        clevs = ticker.MaxNLocator(nbins=256).tick_values(vmin, vmax)
        outfile = os.path.join(out_path,'%s_%s_%s.png' %(region, variable, trec.strftime('%Y%m%d%H')))
        plot_tri_var(grid, val, trec, plotbound, clevs, 'rainbow', vmin, vmax, stride, 
                     title, units, cb_location, outfile, region)
    


