#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:46:14 2022

@author: sg
"""
import numpy as np
import inputs as inp
import datamap_plotter as plotter
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import utils

gname = '2A0335'
telescope = 'MUSE'
# telescope = 'ALMA'
# pltmap = 'velocity'
pltmap = 'flux'
max_error = 5

sysparam = inp.get_sys_params(gname,telescope,pltmap)

if telescope=='MUSE':
    if pltmap=='velocity':
        fnvel = sysparam['vmapfn']
        fnvelerr = sysparam['verrmapfn']
        fnflux = sysparam['fluxmapfn']
        flcut = sysparam['flcut']
        ra = sysparam['ra']
        dec = sysparam['dec']
        res = sysparam['res']
        nx, ny = sysparam['nx'], sysparam['ny']
        ylab = sysparam['ylab']
        xlab = sysparam['xlab']
        msk_sz = sysparam['msk_sz']
        rkpc = sysparam['rkpc']
        xpl, xpu = sysparam['xpl'], sysparam['xpu']
        ypl, ypu = sysparam['ypl'], sysparam['ypu']
        cbfrac = sysparam['cbfrac']
        one_third, half = sysparam['one_third'], sysparam['half']
        vsfyl, vsfyu = sysparam['vsfyl'], sysparam['vsfyu']
        histulim, sepbin = sysparam['histulim'], sysparam['sepbin']
        res = res*0.2  # constant sampling size is 0.2"
    if pltmap=='flux':
        fnflux = sysparam['fluxmapfn']
        flcut = sysparam['flcut']
        ra = sysparam['ra']
        dec = sysparam['dec']
        res = sysparam['res']
        nx, ny = sysparam['nx'], sysparam['ny']
        ylab = sysparam['ylab']
        xlab = sysparam['xlab']
        msk_sz = sysparam['msk_sz']
        rkpc = sysparam['rkpc']
        xpl, xpu = sysparam['xpl'], sysparam['xpu']
        ypl, ypu = sysparam['ypl'], sysparam['ypu']
        cbfrac = sysparam['cbfrac']
        one_third, half = sysparam['one_third'], sysparam['half']
        vsfyl, vsfyu = sysparam['vsfyl'], sysparam['vsfyu']
        histulim, sepbin = sysparam['histulim'], sysparam['sepbin']
        res = res*0.2  # constant sampling size is 0.2"
if telescope=='ALMA':
    if pltmap=='velocity':
        fnvel = sysparam['vmapfn']
        ra = sysparam['ra']
        dec = sysparam['dec']
        res = sysparam['res']
        nx, ny = sysparam['nx'], sysparam['ny']
        ylab = sysparam['ylab']
        xlab = sysparam['xlab']
        msk_sz = sysparam['msk_sz']
        rkpc = sysparam['rkpc']
        xpl, xpu = sysparam['xpl'], sysparam['xpu']
        ypl, ypu = sysparam['ypl'], sysparam['ypu']
        cbfrac = sysparam['cbfrac']
        one_third, half = sysparam['one_third'], sysparam['half']
        vsfyl, vsfyu = sysparam['vsfyl'], sysparam['vsfyu']
        histulim, sepbin = sysparam['histulim'], sysparam['sepbin']
        n_bins, cd = sysparam['n_bins'], sysparam['cd']
        res = res*cd # sampling size changes for each source
    # if pltmap=='flux':

## extracting BH location
if pltmap=='velocity':
    filename = get_pkg_data_filename(fnvel)
if pltmap=='flux':
    filename = get_pkg_data_filename(fnflux)
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
c = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.deg))
target = utils.skycoord_to_pixel(c, wcs)
Cx, Cy = target[0],target[1]
print(Cx,Cy)

## plotting the velocity maps with/without masks
if telescope=='MUSE':
    if pltmap=='velocity':
        velo = fits.open(fnvel)
        velo_data = velo[0].data
        velo.close()
        error = fits.open(fnvelerr)
        error_data = error[0].data
        error.close()
        flux = fits.open(fnflux)
        flux_data = flux[0].data
        flux.close()
        alldata = {'gname':gname,'velo_data':velo_data,'error_data':error_data,
                   'flux_data':flux_data,'flcut':flcut,'rkpc':rkpc,'mask_size':msk_sz,
                   'Cx':Cx,'Cy':Cy,'ylab':ylab,'xlab':xlab,'res':res,
                   'xpl':xpl, 'xpu':xpu, 'ypl':ypl, 'ypu':ypu, 'cbfrac':cbfrac}
        good_v = plotter.plot_data(alldata, max_error, telescope, pltmap)
    if pltmap=='flux':
        flux = fits.open(fnflux)
        flux_data = flux[0].data
        flux.close()
        alldata = {'gname':gname,'flux_data':flux_data,'flcut':flcut,
                   'mask_size':msk_sz,'Cx':Cx,'Cy':Cy,'ylab':ylab,'xlab':xlab,
                   'res':res,'xpl':xpl,'xpu':xpu,'ypl':ypl,'ypu':ypu,'cbfrac':cbfrac}
        good_v = plotter.plot_data(alldata, max_error, telescope, pltmap)