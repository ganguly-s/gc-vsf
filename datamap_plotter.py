#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 07:37:04 2022

@author: sg
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import datamap_mask as mask

#def plot_data(alldata, max_error, telescope, pltmap):
def plot_velo_data(all_data, max_error, vrange, cuts, flux_cut=True, rand_mask=True):
    # extracting plot specifications
    Cx, Cy = all_data['Cx'], all_data['Cy']
    rkpc, res = all_data['rkpc'], all_data['res']
    gname = all_data['gname']
    ylab, xlab = all_data['ylab'], all_data['xlab']
    xpl, xpu, ypl, ypu, cbfrac = all_data['xpl'], all_data['xpu'], all_data['ypl'], all_data['ypu'], all_data['cbfrac']
    # saving the images
    path0 = 'MUSE/velocity/'
    path = path0+gname+'/'
    if not os.path.exists('MUSE/'):
        os.mkdir('MUSE/')
        if not os.path.exists(path0):
            os.mkdir(path0)
            if not os.path.exists(path):
                os.mkdir(path)
    maskmap = path+gname+'_vdatuse.png'
    velomap = path+gname+'_vmap.png'
    # obtaining masked data for plotting and comparison
    velo_plot, usedata = mask.apply_mask(all_data, max_error, vrange, cuts, flux_cut, rand_mask)
    print("total number of points used:", usedata)
    # velocity map with and without any mask
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize = (16,7)) 
    for axs in ax:
        axs.plot(Cx, Cy, marker="X", color="k", markersize=10,
              linestyle="None", label="Black Hole Position")
        r = rkpc/res
        circle= plt.Circle((Cx, Cy), radius= r, color='k', fill=False)
        axs.add_artist(circle)
    vminimum = np.nanmin(velo_plot)
    vmaximum = np.nanmax(velo_plot) 
    im1 = ax[0].imshow(all_data['velo_data'], vmin=vminimum, vmax=vmaximum, cmap="Spectral_r")
    cb1=plt.colorbar(im1, ax=ax[0],fraction=cbfrac, pad=0.04)
    cb1.ax.tick_params(axis='y', direction='in')
    im2 = ax[1].imshow(velo_plot, cmap="Spectral_r")
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
    totdata = len(all_data['velo_data'][~np.isnan(all_data['velo_data'])])
    print(totdata,usedata)
    ax[0].set_title('all data = %d'%totdata, size=20)
    if flux_cut and rand_mask:
        ax[1].set_title('random mask (%d)+flux cut (%d) = %d'%(cuts[0],cuts[1],usedata), size=20)
    if not(flux_cut) and rand_mask:
        ax[1].set_title('random mask (%d) = %d'%(cuts,usedata), size=20)
    if flux_cut and not(rand_mask):
        ax[1].set_title('flux cut (%d) = %d'%(cuts,usedata), size=20)
    fig.suptitle(gname, size=24)
    plt.tight_layout()
    plt.savefig(maskmap)
    return

#def plot_flux_data():
#def vsf_plotter():