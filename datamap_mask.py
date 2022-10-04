#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:41:06 2022

@author: sg
"""
import numpy as np

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
    masked_v = np.ma.masked_array(velo_data, mask=~good_v)
    return masked_v, len(velo_data[good_v])