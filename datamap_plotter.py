#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 07:37:04 2022

@author: sg
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import datamap_mask as mask
import inputs as inp

newparams = {'font.family':'serif', 'axes.labelsize': 22, 'axes.linewidth': 1, 
              'lines.linewidth': 2, 'figure.figsize': (10, 8),
              'figure.subplot.wspace': 0.,
              'ytick.labelsize': 22, 'xtick.labelsize': 22,
              'ytick.major.pad': 5, 'xtick.major.pad': 5,
              'legend.fontsize': 20, 'legend.frameon': True,
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
        axs.plot(Cx, Cy, marker="X", color="k", markersize=10,linestyle="None")
        r = rkpc/res
        circle= plt.Circle((Cx, Cy), radius= r, color='k', ls='--', fill=False)
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
    circ = plt.scatter([],[], ls='--',edgecolor='k',facecolor='none',label='r = %d kpc'%rkpc)
    ax[0].legend(handles=[circ],loc="upper left", prop={'size': 22},frameon=False)
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
    if not(flux_cut) and not(rand_mask):
        ax[1].set_title('error cut only = %d'%usedata,size=20)
    fig.suptitle(gname, size=24)
    plt.tight_layout()
    plt.savefig(maskmap)
    plt.show()
    
    # single velocity map for the paper
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize = (10,8)) 
    ax.plot(Cx, Cy, marker="X", color="k", markersize=10,linestyle="None")
    r = rkpc/res
    circle= plt.Circle((Cx, Cy), radius= r, color='k', ls='--', fill=False)
    ax.add_artist(circle)
    vminimum = np.nanmin(velo_plot)
    vmaximum = np.nanmax(velo_plot) 
    im = ax.imshow(velo_plot, cmap="Spectral_r")
    cb=plt.colorbar(im, ax=ax,fraction=cbfrac, pad=0.04)
    cb.ax.tick_params(axis='y', direction='in')
    cb.set_label("line-of-sight velocity (km/s)")
    ax.set_xlim(xpl/res,xpu/res)
    ax.set_ylim(ypu/res,ypl/res)
    ax.set_xticks(x_locs)
    ax.set_yticks(y_locs)
    ax.set_xticklabels(x_labels)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("x (kpc)")
    ax.invert_yaxis()
    ax.tick_params(which='both',direction='in')
    circ = plt.scatter([],[], ls='--',edgecolor='k',facecolor='none',label='r = %d kpc'%rkpc)
    ax.legend(handles=[circ],loc="upper left", prop={'size': 22},frameon=False)
    ax.set_ylabel("y (kpc)")
    plt.suptitle(gname, size=24)
    plt.tight_layout()
    plt.savefig(velomap)
    plt.show()

    return

def velo_data_panel(all_data, max_error, vrange, cuts, flux_cut=True, rand_mask=True):
    # extracting plot specifications
    Cx, Cy = all_data['Cx'], all_data['Cy']
    rkpc, res = all_data['rkpc'], all_data['res']
    gname = all_data['gname']
    ylab, xlab = all_data['ylab'], all_data['xlab']
    xpl, xpu, ypl, ypu, cbfrac, hcbfrac = all_data['xpl'], all_data['xpu'], all_data['ypl'], all_data['ypu'], all_data['cbfrac'], all_data['hcbfrac']
    # saving the images
    path0 = 'MUSE/velocity/'
    path = path0+gname+'/'
    if not os.path.exists('MUSE/'):
        os.mkdir('MUSE/')
    if not os.path.exists(path0):
        os.mkdir(path0)
    if not os.path.exists(path):
        os.mkdir(path)
    maskmap = path+gname+'_vpanel.pdf'
    # obtaining masked data for plotting
    good_v = mask.apply_mask(all_data, max_error, vrange, cuts, flux_cut, rand_mask)
    velo_plot = np.ma.masked_array(all_data['velo_data'], mask=~good_v)
    print("total number of points used:", len(all_data['velo_data'][good_v]))
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize = (10,8))
    plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
    ax.plot(Cx, Cy, marker="X", color="k", markersize=10,linestyle="None")
    r = rkpc/res
    circle= plt.Circle((Cx, Cy), radius= r, color='k', ls='--', fill=False)
    ax.add_artist(circle)
    vminimum = np.nanmin(velo_plot)
    vmaximum = np.nanmax(velo_plot) 
    im = ax.imshow(velo_plot, cmap="Spectral_r")
    # calculate width/height of image
    cb=plt.colorbar(im, ax=ax,fraction=hcbfrac, pad=0.05,orientation='horizontal')
    cb.ax.tick_params(axis='x', direction='in')
    cb.set_label(r"(km s$^{-1}$)",labelpad=-50)
    ax.set_xlim(xpl/res,xpu/res)
    ax.set_ylim(ypu/res,ypl/res)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.tick_params(which='both',direction='in')
    # adding custom scale bar
    scale = 10#10/res/2
    scaleres = 2*scale # pixel size for MUSE in arcsec and kpc below, factor 2 due to errorbar
    plt.errorbar( 250, 20, xerr=scale, color='k', capsize=5)
    #plt.text( 250, 28, r'%.1f $^{\prime\prime}$'%scaleres,  horizontalalignment='center', verticalalignment='top')
    plt.text( 250, 15, r'%.1f kpc'%(scaleres*res),  horizontalalignment='center', verticalalignment='top',size=16)
    plt.text(20,210,gname,size=20)
    #plt.tight_layout()
    plt.savefig(maskmap,bbox_inches='tight')
    plt.show()

    return

def xray_plotter(all_data,fnflux):
    # extracting plot specifications
    Cx, Cy = all_data['Cx'], all_data['Cy']
    rkpc, res = all_data['rkpc'], all_data['res']
    gname = all_data['gname']
    ylab, xlab = all_data['ylab'], all_data['xlab']
    xpl, xpu, ypl, ypu, cbfrac, hcbfrac = all_data['xpl'], all_data['xpu'], all_data['ypl'], all_data['ypu'], all_data['cbfrac'], all_data['hcbfrac']
    
    resarcs = 0.492
    res = res*resarcs/0.2  # removing sampling size for MUSE and multiplying with that of Chandra
    impath = '../X-ray images/'
    outpath = 'Chandra/'
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outpath += gname+'_xray.pdf'
    fname, xlab, ylab, xpl, xpu, ypl, ypu, cbfrac = inp.get_xray_file(gname)
    fn = impath + fname
    #gc = aplpy.FITSFigure(fn)
    #gc.show_grayscale()
    #gc.show_contour(fnvel,colors='red',levels=1)
    #gc.save(outpath+gname+'.pdf')

    xray = fits.open(fn)
    xray_data = xray[0].data
    xray.close()
    hdu0 = fits.open(fnflux)[0]
    wcs0 = WCS(hdu0.header)
    hdu = fits.open(fn)[0]
    wcs = WCS(hdu.header)
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize = (10,8))
    plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
    ax = plt.subplot(projection=wcs)
    iminimum = np.nanmin(xray_data)
    imaximum = np.nanmax(xray_data)
    ax.plot(Cx, Cy, marker="X", color="k", markersize=5,linestyle="None", transform=ax.get_transform(wcs0))
    r = rkpc/res/0.2*resarcs
    circle= plt.Circle((Cx, Cy), radius= r, color='k', ls='--', fill=False, transform=ax.get_transform(wcs0))
    ax.add_artist(circle)
    im = ax.imshow(xray_data,vmin=iminimum, vmax=2e-6, cmap="Spectral_r")    
    ax.contour(hdu0.data, transform=ax.get_transform(wcs0),levels=[100,500],colors='red',alpha=0.8,linewidths=1)
    lon, lat = ax.coords[0], ax.coords[1]
    ax.set_xlim(xpl,xpu)
    ax.set_ylim(ypu,ypl)
    lon.set_ticks_visible(False)
    lat.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')
    ax.invert_yaxis()
    cb = plt.colorbar(im, ax=ax,fraction=cbfrac, pad=0.05, orientation = 'horizontal',ticks=[0.0e-6,0.5e-6,1e-6,1.5e-6,2e-6])
    ax.tick_params(which='both',direction='in')
    cb.ax.tick_params(axis='x', direction='in')
    cb.set_label(r"($10^{-6}$ ergs s$^{-1}$ cm$^{-2}$)",labelpad=-50)
    cb.ax.set_xticklabels(['0.0','0.5','1.0','1.5','2.0'])
    # adding custom scale bar
    scale = 10#10/res/2
    scaleres = 2*scale # for Chandra
    ax.errorbar( 2060, 2470, xerr=scale, color='w', capsize=5)
    #plt.text( 250, 28, r'%.1f $^{\prime\prime}$'%scaleres,  horizontalalignment='center', verticalalignment='top')
    plt.text( 2060, 2460, r'%.1f kpc'%(scaleres*res),  horizontalalignment='center', verticalalignment='top',color='white',size=16)
    #plt.tight_layout()
    plt.savefig(outpath,bbox_inches='tight')
    plt.show()

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
    fnvsfall = impath+gname+'_allvsf.pdf'
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
    plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
    plt.loglog(dist_array_kpc, y_expect, "k--", label="1/3")
    plt.loglog(dist_array_kpc, y_expect2, "m-.", label="1/2")
    plt.loglog(dist_array_kpc, v_diff_mean_smooth, marker="o",linestyle="None", markersize=4, color="C0",label='All')
    plt.fill_between(dist_array_kpc[:-1], lower_v[:-1], upper_v[:-1], color="C0", alpha=0.2)
    plt.loglog(dist_array_kpc_in, v_diff_mean_smooth_in, marker="s",linestyle="None", markersize=4, color="C3",label='inner (r<'+str(rkpc)+' kpc)')
    plt.fill_between(dist_array_kpc_in[:-1], lower_v_in[:-1], upper_v_in[:-1], color="C3", alpha=0.2)
    plt.loglog(dist_array_kpc_out, v_diff_mean_smooth_out, marker="D",linestyle="None", markersize=4, color="C2",label='outer (r>'+str(rkpc)+' kpc)')
    plt.fill_between(dist_array_kpc_out[:-1], lower_v_out[:-1], upper_v_out[:-1], color="C2", alpha=0.2)
    plt.xlabel(r"$\mathcal{l}$ (kpc)",size=22)
    plt.ylabel(r"$<|\delta v|>\, \rm (km/s)$",size=22)
    plt.legend(loc="lower right",fontsize=20)
    plt.ylim(vsfyl, vsfyu)
    plt.grid()
    plt.tick_params(which='both',direction='in')
    #plt.title(gname, size=24)
    plt.savefig(fnvsfall,bbox_inches='tight')

def bplplotter(telescope,gname,vsfyl,vsfyu,one_third,half,pars1,stdevs1,pars2,stdevs2,b1,b2):
    impath = telescope+'/vsffits/'+gname+'/'
    if not os.path.exists(telescope+'/'):
        os.mkdir(telescope+'/')
    if not os.path.exists(telescope+'/vsffits/'):
        os.mkdir(telescope+'/vsffits/')
    if not os.path.exists(telescope+'/vsffits/'+gname+'/'):
        os.mkdir(telescope+'/vsffits/'+gname+'/')
    datpath = telescope+'/vsfdat/'+gname+'/'
    fnvsfall = impath+gname+'_allvsf.png'
    fnvsfallnpz = datpath+gname+'_allvsf.npz'

    vsfdata= np.load(fnvsfallnpz)
    rkpc = vsfdata['rkpc']
    dist_array = vsfdata['dist_array']
    dist_array_kpc = vsfdata['dist_array_kpc']
    v_diff_mean_smooth = vsfdata['v_diff_mean_smooth']

    y_expect = dist_array**(1.0/3)*one_third
    y_expect2 = dist_array**(1.0/2)*half
    
    plt.figure(figsize=(10, 8))
    plt.loglog(dist_array_kpc, y_expect, "k--", label="1/3")
    plt.loglog(dist_array_kpc, y_expect2, "m-.", label="1/2")
    plt.loglog(dist_array_kpc, v_diff_mean_smooth, marker="o",linestyle="None", markersize=4, color="C0",label=gname+r' All, $b_1$ = %.2f'%b1)

    for i in range(len(dist_array_kpc)):
        if dist_array_kpc[i]>b1:
            distspec = dist_array_kpc[0:i+1]
            break
    plt.loglog(distspec, 10**(pars1[0]+pars1[1]*np.log10(distspec)), 'x', color='C3',markersize=4, label=r'$\alpha_1$ = %.2f $\pm$ %.3f'%(pars1[1],stdevs1[1]))
    if b2!=0:
        for j in range(i,len(dist_array_kpc)):
            if dist_array_kpc[j]>b2:
                distspec = dist_array_kpc[i:j+1]
                break
        plt.loglog(distspec, 10**(pars2[0]+pars2[1]*np.log10(distspec)), 'x', color='C2',markersize=4, label=r'$\alpha_2$ = %.2f $\pm$ %.3f'%(pars2[1],stdevs2[1]))
    plt.xlabel("separation (kpc)")
    plt.ylabel(r"$|\delta v|\, \rm (km/s)$")
    plt.legend(loc="lower right", prop={'size': 20})
    plt.ylim(vsfyl, vsfyu)
    plt.grid()
    plt.tick_params(which='both',direction='in')
    plt.title(gname, size=24)
    plt.savefig(fnvsfall)
    plt.show()
    
    return
