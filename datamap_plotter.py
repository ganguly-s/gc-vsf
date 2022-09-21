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
    
def make_data(alldata,good_v,pltmap):
    if pltmap=='velocity':
        data, error_data = alldata['velo_plot'], alldata['error_data']
        Cx, Cy = alldata['Cx'], alldata['Cy']
        rkpc, res = alldata['rkpc'], alldata['res']
        r = rkpc/res
    # begin calculations of l and delta_v
    xvalues = np.arange(0, len(data[0, :]))
    yvalues = np.arange(0, len(data[:, 0]))
    xx, yy = np.meshgrid(xvalues, yvalues)
    
    goodxx = xx[good_v]
    goodyy = yy[good_v]

    #good_velo = make_mask(velo_data[good_v])
    vel_a = np.reshape(data[good_v], (data[good_v].size, 1))
    vel_b = np.reshape(data[good_v], (1, data[good_v].size))
    v_diff_matrix = np.subtract(
        vel_a, vel_b, dtype=np.float64)  # vel_a - vel_b

    error_a = np.reshape(error_data[good_v]**2, (error_data[good_v].size, 1))
    error_b = np.reshape(error_data[good_v]**2, (1, error_data[good_v].size))
    error_matrix = error_a + error_b

    px_a = np.reshape(goodxx, (goodxx.size, 1))
    px_b = np.reshape(goodxx, (1, goodxx.size))
    py_a = np.reshape(goodyy, (goodyy.size, 1))
    py_b = np.reshape(goodyy, (1, goodyy.size))
    dist_matrix = np.sqrt((px_a - px_b)**2 + (py_a - py_b)**2)
    v_diff_half = np.ndarray.flatten(np.triu(v_diff_matrix, k=0))
    # this is still a 2D matrix, with the lower half values all set to 0
    dist_half = np.ndarray.flatten(np.triu(dist_matrix, k=0))
    error_half = np.ndarray.flatten(np.triu(error_matrix, k=0))
    good_dist = dist_half > 0
    np_dist = dist_half[good_dist]
    np_v_diff = v_diff_half[good_dist]
    np_error_2 = error_half[good_dist]
# np.savez("M87_save.npz",np_dist=np_dist,np_v_diff=np_v_diff,np_error_2=np_error_2)
    acx=np.zeros((1,goodxx.size))+Cx
    acy=np.zeros((1,goodxx.size))+Cy
    bcx=np.zeros((goodxx.size,1))+Cx
    bcy=np.zeros((goodxx.size,1))+Cy

    rA_matrix=np.sqrt((px_a - acx)**2 + (py_a - acy)**2)
    rB_matrix=np.sqrt((px_b - bcx)**2 + (py_b - bcy)**2)
    rA_half = np.ndarray.flatten(np.triu(rA_matrix, k=0)) # this is still a 2D matrix, with the lower half values all set to 0 
    rB_half = np.ndarray.flatten(np.triu(rB_matrix, k=0)) # this is still a 2D matrix, with the lower half values all set to 0 
    rA=rA_half[good_dist]
    rB=rB_half[good_dist]
    inner_f=(rA<=r)&(rB<=r)
    outer_f=(rA>r)&(rB>r)
    np_dist_in=np_dist[inner_f]
    np_v_diff_in=np_v_diff[inner_f]
    np_error_in_2=np_error_2[inner_f]
    np_dist_out=np_dist[outer_f]
    np_v_diff_out=np.abs(np_v_diff[outer_f])
    # np_v_diff_2=np_v_diff**2
    # np_v_diff_3=np_v_diff**3
    # np_v_diff_4=np_v_diff**4
    np_error_out_2=np_error_2[outer_f]
    return np_dist, np_v_diff, np_error_2, np_dist_in, np_v_diff_in, np_error_in_2, np_dist_out, np_v_diff_out, np_error_out_2