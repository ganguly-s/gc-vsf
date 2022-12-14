#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:41:06 2022

@author: sg
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# create a mask to select a fraction of the points.
def make_mask(velo_data, mask_size):
    a = np.zeros(velo_data.shape, dtype=int)
    b = a.flatten()
    b[0:mask_size] = 1
    np.random.shuffle(b)
    b = b.astype(bool)
    a = b.reshape(velo_data.shape)
    return a

# removes nans for smoothing final VSF
def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]

def run_helper(x, y):
    nans, x = nan_helper(y)
    z = np.copy(y)
    z[nans] = np.interp(x(nans), x(~nans), y[~nans])
    return(z)

# apply masks according to user choice
def apply_mask(all_data, max_error, vrange, cuts, flux_cut, rand_mask):
    velo_data = all_data['velo_data']
    error_data = all_data['error_data']
    flux_data = all_data['flux_data']
    if flux_cut and rand_mask:
        msk_size, flcut = cuts[0], cuts[1]
        this_mask=make_mask(velo_data=velo_data,mask_size=msk_size)
        good_v = (velo_data < vrange) & (velo_data > -vrange) & (error_data < max_error) & (flux_data > flcut) & this_mask
    if not(flux_cut) and rand_mask:
        msk_size = cuts
        this_mask=make_mask(velo_data=velo_data,mask_size=msk_size)
        good_v = (velo_data < vrange) & (velo_data > -vrange) & (error_data < max_error) & this_mask
    if flux_cut and not(rand_mask):
        flcut = cuts
        good_v = (velo_data < vrange) & (velo_data > -vrange) & (error_data < max_error) & (flux_data > flcut)
    if not(flux_cut) and not(rand_mask):
        good_v = (velo_data < vrange) & (velo_data > -vrange) & (error_data < max_error)
    
    return good_v

def power_law(x,a,b):
    return a+(x*b)

def brokenpowerlaw(telescope,gname,r1,r2):
    #impath = telescope+'/vsfplots/'+gname+'/'
    datpath = telescope+'/vsfdat/'+gname+'/'
    #fnvsfall = impath+gname+'_allvsf.png'
    fnvsfallnpz = datpath+gname+'_allvsf.npz'
    
    vsfdata= np.load(fnvsfallnpz)
    rkpc = vsfdata['rkpc']
    dist_array = vsfdata['dist_array']
    dist = vsfdata['dist_array_kpc']
    vsf = vsfdata['v_diff_mean_smooth']
    
    for i in range(len(dist)):
        if dist[i]>r1:
            b1 = i
            break
    
    # fitting upto the first break
    distspec = dist[0:b1+1]
    vsfspec = vsf[0:b1+1]
    
    pars1, cov1 = curve_fit(f=power_law, xdata=np.log10(distspec), ydata=np.log10(vsfspec), p0=[10., 0.33])
    stdevs1 = np.sqrt(np.diag(cov1))

    # fitting from first break upto the viable peak
    if r2!=0: # for single power law, it is
        for i in range(b1,len(dist)):
            if dist[i]>r2:
                b2 = i
                break
        distspec = dist[b1:b2+1]
        vsfspec = vsf[b1:b2+1]

        pars2, cov2 = curve_fit(f=power_law, xdata=np.log10(distspec), ydata=np.log10(vsfspec), p0=[10., 0.33])
        stdevs2 = np.sqrt(np.diag(cov2))
    else:
        pars2,stdevs2 = [], []

    return pars1,stdevs1,pars2,stdevs2
