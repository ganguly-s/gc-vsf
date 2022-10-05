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

newparams = {'font.family':'serif', 'axes.labelsize': 24, 'axes.linewidth': 1, 
              'lines.linewidth': 2, 'figure.figsize': (8, 6),
              'figure.subplot.wspace': 0.4,
              'ytick.labelsize': 24, 'xtick.labelsize': 24,
              'ytick.major.pad': 5, 'xtick.major.pad': 5,
              'legend.fontsize': 24, 'legend.frameon': True,
              'legend.handlelength': 2, 'figure.facecolor':'white',
              'savefig.dpi': 300}
plt.rcParams.update(newparams)

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
    good_v = mask.apply_mask(all_data, max_error, vrange, cuts, flux_cut, rand_mask)
    velo_plot = np.ma.masked_array(all_data['velo_data'], mask=~good_v)
    print("total number of points used:", len(all_data['velo_data'][good_v]))
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
    y_labels = ylab
    y_locs = y_labels/res
    x_labels = xlab
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
    usedata = len(all_data['velo_data'][good_v])
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
def vsf_plotter(gname,telescope,vsfyl,vsfyu,one_third,half):
    impath = telescope+'/vsfplots/'+gname+'/'
    datpath = telescope+'/vsfdat/'+gname+'/'
    if not os.path.exists(telescope+'/'):
        os.mkdir(telescope+'/')
    if not os.path.exists(telescope+'/vsfplots/'):
        os.mkdir(telescope+'/vsfplots/')
    if not os.path.exists(telescope+'/vsfplots/'+gname+'/'):
        os.mkdir(telescope+'/vsfplots/'+gname+'/')
    fnvsfall = impath+gname+'_allvsf.png'
    fnvsfallnpz = datpath+gname+'_allvsf.npz'
    
    vsfdata= np.load(fnvsfallnpz)
    rkpc = vsfdata['rkpc']
    dist_array = vsfdata['dist_array']
    dist_array_kpc = vsfdata['dist_array_kpc']
    dist_array_kpc_in = vsfdata['dist_array_kpc_in']
    dist_array_kpc_out = vsfdata['dist_array_kpc_out']
    v_diff_mean_smooth = vsfdata['v_diff_mean_smooth']
    v_diff_mean_smooth_in = vsfdata['v_diff_mean_smooth_in']
    v_diff_mean_smooth_out = vsfdata['v_diff_mean_smooth_out']
    lower_v = vsfdata['lower_v']
    lower_v_in = vsfdata['lower_v_in']
    lower_v_out = vsfdata['lower_v_out']
    upper_v = vsfdata['upper_v']
    upper_v_in = vsfdata['upper_v_in']
    upper_v_out = vsfdata['upper_v_out']
    
    y_expect = dist_array**(1.0/3)*one_third
    y_expect2 = dist_array**(1.0/2)*half
    
    plt.figure(figsize=(10, 8))
    plt.loglog(dist_array_kpc, y_expect, "k--", label="1/3")
    plt.loglog(dist_array_kpc, y_expect2, "m-.", label="1/2")
    plt.loglog(dist_array_kpc, v_diff_mean_smooth, marker="o",linestyle="None", markersize=4, color="C0",label=gname+' All')
    plt.fill_between(dist_array_kpc[:-1], lower_v[:-1], upper_v[:-1], color="C0", alpha=0.2)
    plt.loglog(dist_array_kpc_in, v_diff_mean_smooth_in, marker="s",linestyle="None", markersize=4, color="C3",label='inner (r<'+str(rkpc)+' kpc)')
    plt.fill_between(dist_array_kpc_in[:-1], lower_v_in[:-1], upper_v_in[:-1], color="C3", alpha=0.2)
    plt.loglog(dist_array_kpc_out, v_diff_mean_smooth_out, marker="D",linestyle="None", markersize=4, color="C2",label='outer (r>'+str(rkpc)+' kpc)')
    plt.fill_between(dist_array_kpc_out[:-1], lower_v_out[:-1], upper_v_out[:-1], color="C2", alpha=0.2)
    plt.xlabel("separation (kpc)")
    plt.ylabel(r"$|\delta v|\, \rm (km/s)$")
    plt.legend(loc="lower right", prop={'size': 22})
    plt.ylim(vsfyl, vsfyu)
    plt.grid()
    plt.tick_params(which='both',direction='in')
    plt.title(gname, size=24)
    plt.savefig(fnvsfall)