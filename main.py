#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:46:14 2022

@author: sg
"""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import utils
import inputs as inp
import datamap_plotter as plotter
import datamap_mask as mask
import datamap_cleaner as clean

gname = '2A0335'
telescope = 'MUSE'
# telescope = 'ALMA'
pltmap = 'velocity'
# pltmap = 'flux'
max_error = 5

sysparam = inp.get_sys_params(gname,telescope,pltmap)

# extracting BH location
def gal_center(fn):
    hdu = fits.open(fn)[0]
    wcs = WCS(hdu.header)
    c = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.deg))
    target = utils.skycoord_to_pixel(c, wcs)
    Cx, Cy = target[0],target[1]
    print("center:",Cx,Cy)
    return Cx, Cy

if telescope=='MUSE':
    if pltmap=='velocity':
        fnvel = 'inputs/'+sysparam['vmapfn']
        fnvelerr = 'inputs/'+sysparam['verrmapfn']
        fnflux = 'inputs/'+sysparam['fluxmapfn']
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
        vsfyl, vsfyu, vsfbin = sysparam['vsfyl'], sysparam['vsfyu'], sysparam['vsfbin']
        histulim, sepbin = sysparam['histulim'], sysparam['sepbin']
        res = res*0.2  # constant sampling size is 0.2"
        # data extraction and plotting
        velo = fits.open(fnvel)
        velo_data = velo[0].data
        velo.close()
        error = fits.open(fnvelerr)
        error_data = error[0].data
        error.close()
        flux = fits.open(fnflux)
        flux_data = flux[0].data
        flux.close()
        Cx, Cy = gal_center(fnvel)
        # setting up for plots
        flux_cut = True
        rand_mask = False
        if flux_cut and rand_mask:
            cuts = [msk_sz, flcut]
        if flux_cut and not(rand_mask):
            cuts = flcut
        if not(flux_cut) and rand_mask:
            cuts = msk_sz
        # Choosing the type of mask to be applied to data, error cut is mandatory, although it barely has an effect
        alldata = {'gname':gname,'velo_data':velo_data, 'error_data':error_data, 'flux_data':flux_data,
                   'Cx':Cx,'Cy':Cy,'rkpc':rkpc,'res':res,'ylab':ylab,'xlab':xlab,
                    'xpl':xpl, 'xpu':xpu, 'ypl':ypl, 'ypu':ypu, 'cbfrac':cbfrac}
        plotter.plot_velo_data(alldata, max_error, 2000, cuts, flux_cut, rand_mask)
        # good_v = mask.apply_mask(alldata, max_error, 2000, cuts, flux_cut, rand_mask)
        # alldata['good_v'] = good_v
        # alldata['telescope'] = telescope
        # alldata['one_third'], alldata['half'] =  one_third, half
        # alldata['vsfyl'], alldata['vsfyu'], alldata['vsfbin'] = vsfyl, vsfyu, vsfbin
        # clean.calc_vsf(alldata)
        plotter.vsf_plotter(gname, telescope, vsfyl, vsfyu, one_third, half)
    if pltmap=='flux':
        fnflux = 'inputs/'+sysparam['fluxmapfn']
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
        # data extraction and plotting
        flux = fits.open(fnflux)
        flux_data = flux[0].data
        flux.close()
        alldata = {'gname':gname,'flux_data':flux_data,'flcut':flcut,
                   'mask_size':msk_sz,'Cx':Cx,'Cy':Cy,'ylab':ylab,'xlab':xlab,
                   'res':res,'xpl':xpl,'xpu':xpu,'ypl':ypl,'ypu':ypu,'cbfrac':cbfrac}
        #plotter.plot_data(alldata, max_error, telescope, pltmap)
if telescope=='ALMA':
    if pltmap=='velocity':
        fnvel = 'inputs/'+sysparam['vmapfn']
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
        