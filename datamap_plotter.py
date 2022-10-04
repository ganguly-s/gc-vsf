#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 07:37:04 2022

@author: sg
"""
import numpy as np
import matplotlib.pyplot as plt

# create a mask to select a fraction of the points.
def make_mask(velo_data, mask_size):
    a = np.zeros(velo_data.shape, dtype=int)
    b = a.flatten()
    b[0:mask_size] = 1
    np.random.shuffle(b)
    b = b.astype(bool)
    a = b.reshape(velo_data.shape)
    return a

def plot_data(alldata, max_error, telescope, pltmap):
    gname = alldata['gname']
    if telescope=='MUSE':
        impath = 'images/'
        if pltmap=='velocity':
            data, error_data = alldata['velo_data'], alldata['error_data']
            fldata, flcut = alldata['flux_data'], alldata['flcut']
            datrangelow, datrangeup = -200, 400
            rkpc = alldata['rkpc']
            maskmap = impath+gname+'_Halpha_vmap.png'
        if pltmap=='flux':
            data = alldata['flux_data']
            datrangelow, datrangeup = 0.001, 1000
            maskmap = impath+gname+'_Halpha_fluxmap.png'
    mask_size, flcut = alldata['mask_size'], alldata['flcut']
    Cx, Cy = alldata['Cx'], alldata['Cy']
    ylab, xlab, res = alldata['ylab'], alldata['xlab'], alldata['res']
    xpl, xpu, ypl, ypu, cbfrac = alldata['xpl'], alldata['xpu'], alldata['ypl'], alldata['ypu'], alldata['cbfrac']
    
    # choosing random numbers from data set
    this_mask=make_mask(data,mask_size)
    
    # setting a range for the data, error and flux cut if available, 
    # in addition to a random data point selection mask
    if pltmap=='velocity':
        good_v = (data < datrangeup) & (data > datrangelow) & (error_data < max_error) & (fldata > flcut) & this_mask
    if pltmap=='flux':
        good_v = (data < np.nanmax(data)) & (data > datrangelow) & (data > flcut) & this_mask
    print("good_v", good_v.shape)
    print("total number of points used:", data[good_v].shape)
    data_plot = np.ma.masked_array(data, mask=~good_v)
    #print(velo_plot.size)
    
    # velocity map with and without any mask
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize = (16,7)) 
    for axs in ax:
        axs.plot(Cx, Cy, marker="X", color="k", markersize=10,
              linestyle="None", label="Black Hole Position")
        if pltmap=='velocity':
            r = rkpc/res
            circle= plt.Circle((Cx, Cy), radius= r, color='k', fill=False)
            axs.add_artist(circle)
    vminimum = np.nanmin(data)
    vmaximum = np.nanmax(data) 
    im1 = ax[0].imshow(data, vmin=vminimum, vmax=vmaximum, cmap="Spectral_r")
    cb1=plt.colorbar(im1, ax=ax[0],fraction=cbfrac, pad=0.04)
    cb1.ax.tick_params(axis='y', direction='in')
    im2 = ax[1].imshow(data_plot, cmap="Spectral_r")
    cb2=plt.colorbar(im2, ax=ax[1],fraction=cbfrac, pad=0.04)
    cb2.ax.tick_params(axis='y', direction='in')
    cb2.set_label("line-of-sight velocity (km/s)")
    y_labels = ylab #np.arange(1, 6)
    y_locs = y_labels/res
    x_labels = xlab #np.arange(0, 6)
    x_locs = x_labels/res
    for axs in ax:
        axs.set_xlim(xpl/res,xpu/res)
        axs.set_ylim(ypu/res,ypl/res)
        axs.set_xticks(x_locs)
        axs.set_yticks(y_locs)
        axs.set_xticklabels(x_labels)
        axs.set_yticklabels(y_labels)
        axs.set_xlabel("x (kpc)")
        axs.invert_yaxis()
        axs.tick_params(which='both',direction='in')
    ax[0].legend(loc="upper left", prop={'size': 22})
    ax[0].set_ylabel("y (kpc)")
    totdata = len(data[~np.isnan(data)])
    usedata = len(data[good_v])
    print(totdata,usedata)
    ax[0].set_title('all data = %d'%totdata, size=20)
    ax[1].set_title('this_mask(%d)+flux_cut(%d) = %d'%(mask_size,flcut,usedata), size=20)
    fig.suptitle(gname, size=24)
    plt.tight_layout()
    plt.savefig(maskmap)
    
    return good_v