#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 08:29:50 2022

@author: sg
"""
import numpy as np
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

def get_sys_params(gname, telescope, pltmap):
    if telescope=='MUSE':
        if gname=='Centaurus':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'Centaurus_velo_s.fits', 'verrmapfn':'Centaurus_velo_err_s.fits',
                           'ra':'12:48:49.16', 'dec':'-41:18:40.9', 'res':0.208, 'nx':359, 'ny':362,
                           'ylab':np.array([5,10,15]), 'xlab':np.array([0,5,10,15]), 'msk_sz':50000,
                           'rkpc':2, 'xpl':0, 'xpu':15, 'ypl':0, 'ypu':15, 'cbfrac':0.040, 'one_third':42, 
                           'half':28, 'vsfyl':7, 'vsfyu':200, 'histulim':80, 'sepbin':300}
            # if pltmap=='flux':
        if gname=='Hydra-A':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'Hydra-A_velo_s.fits', 'verrmapfn':'Hydra-A_velo_err_s.fits',
                           'ra':'09:18:05.68', 'dec':'−12:05:43.8', 'res':1.053, 'nx':361, 'ny':362,
                           'ylab':np.array([25,30,40,50]), 'xlab':np.array([20,30,40,50]), 
                           'msk_sz':1000000, 'rkpc':5, 'xpl':20, 'xpu':50, 'ypl':25, 'ypu':55, 
                           'cbfrac':0.040, 'one_third':88, 'half':55, 'vsfyl':10, 'vsfyu':400, 
                           'histulim':80, 'sepbin':300}
            # if pltmap=='flux':
        if gname=='2A0335':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'2A0335_velo.fits', 'verrmapfn':'2A0335_velo_err.fits',
                           'ra':'03:38:40.54612', 'dec':'+09:58:12.1373', 'res':0.700, 'nx':360, 'ny':363,
                           'ylab':np.array([10,20,30]), 'xlab':np.array([0,10,20,30,40]), 'vsfbin':200,
                           'msk_sz':70000, 'rkpc':5, 'xpl':0, 'xpu':40, 'ypl':0, 'ypu':35, 
                           'cbfrac':0.040, 'one_third':42, 'half':25, 'vsfyl':10, 'vsfyu':300, 
                           'histulim':80, 'sepbin':300, 'fluxmapfn':'2A0335_Halpha_flux.fits','flcut':250}
            if pltmap=='flux':
                sys_par = {'fluxmapfn':'2A0335_Halpha_flux.fits', 'flcut':250, 'ra':'03:38:40.54612', 
                           'dec':'+09:58:12.1373', 'res':0.700, 'nx':360, 'ny':363,
                           'ylab':np.array([10,20,30]), 'xlab':np.array([0,10,20,30,40]), 
                           'msk_sz':90000, 'rkpc':5, 'xpl':0, 'xpu':40, 'ypl':0, 'ypu':35, 
                           'cbfrac':0.040, 'one_third':42, 'half':25, 'vsfyl':10, 'vsfyu':300, 
                           'histulim':80, 'sepbin':300}
        if gname=='A1795':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'A1795_velo.fits', 'verrmapfn':'A1795_velo_err.fits',
                           'ra':'13:48:52.495', 'dec':'+26:35:34.32', 'res':1.220, 'nx':360, 'ny':362,
                           'ylab':np.array([0,10,20,30,40,50,60,70,80]), 'xlab':np.array([20,30,40,50,60,70]), 
                           'msk_sz':1000000, 'rkpc':8, 'xpl':20, 'xpu':70, 'ypl':0, 'ypu':80, 
                           'cbfrac':0.040, 'one_third':55, 'half':33, 'vsfyl':10, 'vsfyu':300, 
                           'histulim':80, 'sepbin':200,'fluxmapfn':'A1795_Halpha_flux.fits',
                           'flcut':250,'vsfbin':200}
            # if pltmap=='flux':
        if gname=='A3581':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'A3581_velo_s.fits', 'verrmapfn':'A3581_velo_err_s.fits',
                           'ra':'14:07:29.76313', 'dec':'-27:01:04.5689', 'res':0.435, 'nx':362, 'ny':361,
                           'ylab':np.array([0,10,20,30]), 'xlab':np.array([10,20,30]), 
                           'msk_sz':1000000, 'rkpc':3, 'xpl':0, 'xpu':30, 'ypl':5, 'ypu':35, 
                           'cbfrac':0.040, 'one_third':60, 'half':33, 'vsfyl':6, 'vsfyu':300, 
                           'histulim':80, 'sepbin':200, 'fluxmapfn':'A3581_Halpha_flux.fits',
                           'flcut':150,'vsfbin':200}
            # if pltmap=='flux':
                
        if gname=='PKS0745':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'pks0745_velo_s.fits', 'verrmapfn':'pks0745_velo_err_s.fits',
                           'ra':'07:47:31.321', 'dec':'-19:17:39.97', 'res':1.890, 'nx':360, 'ny':362,
                           'ylab':np.array([60,70,80,90,100]), 'xlab':np.array([55,70,80,90,100,110]), 
                           'msk_sz':1000000, 'rkpc':3, 'xpl':55, 'xpu':110, 'ypl':55, 'ypu':100, 
                           'cbfrac':0.038, 'one_third':35, 'half':25, 'vsfyl':9, 'vsfyu':150, 
                           'histulim':80, 'sepbin':200}
            # if pltmap=='flux':
        if gname=='R0821':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'RXJ0821_velo_s.fits', 'verrmapfn':'RXJ0821_velo_err_s.fits',
                           'ra':'08:21:02.25', 'dec':'+07:51:47.4', 'res':2.000, 'nx':359, 'ny':363,
                           'ylab':np.array([65,75,85]), 'xlab':np.array([70,80,90,100]), 
                           'msk_sz':1000000, 'rkpc':5, 'xpl':65, 'xpu':100, 'ypl':65, 'ypu':90, 
                           'cbfrac':0.033, 'one_third':42, 'half':28, 'vsfyl':7, 'vsfyu':125, 
                           'histulim':80, 'sepbin':100}
            # if pltmap=='flux':
        if gname=='R1539':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'RXJ1539_velo_s.fits', 'verrmapfn':'RXJ1539_velo_err_s.fits',
                           'ra':'15:39:34.47', 'dec':'-83:35:22.1', 'res':1.437, 'nx':361, 'ny':368,
                           'ylab':np.array([20,40,60,80,100]), 'xlab':np.array([0,20,40,60,80,100]), 
                           'msk_sz':10000, 'rkpc':15, 'xpl':0, 'xpu':100, 'ypl':0, 'ypu':110, 
                           'cbfrac':0.033, 'one_third':37, 'half':20, 'vsfyl':7, 'vsfyu':300, 
                           'histulim':80, 'sepbin':100}
            # if pltmap=='flux':
        if gname=='S1101':
            if pltmap=='velocity':
                sys_par = {'vmapfn':'S1101_velo_s.fits', 'verrmapfn':'S1101_velo_err_s.fits',
                           'ra':'23:13:58.65', 'dec':'-42:43:39.6', 'res':1.094, 'nx':362, 'ny':364,
                           'ylab':np.array([30,40,50,60,70]), 'xlab':np.array([30,40,50,60]), 
                           'msk_sz':1000000, 'rkpc':5, 'xpl':25, 'xpu':60, 'ypl':25, 'ypu':80, 
                           'cbfrac':0.040, 'one_third':45, 'half':25, 'vsfyl':7, 'vsfyu':200, 
                           'histulim':80, 'sepbin':200}
            # if pltmap=='flux':
    if telescope=='ALMA':
        if gname=='Centaurus':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('centaurus_30km_mom1_new.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'centaurus_30km_mom1_new.fits', 'ra':'12:48:49.16', 
                           'dec':'-41:18:40.9', 'res':0.208, 'nx':181, 'ny':181,
                           'ylab':np.array([2,4,6,8]), 'xlab':np.array([0,2,4,6,8]), 
                           'msk_sz':50000, 'rkpc':2, 'xpl':0, 'xpu':8, 'ypl':0, 'ypu':8, 
                           'cbfrac':0.040, 'one_third':62, 'half':48, 'vsfyl':7, 'vsfyu':300, 
                           'histulim':80, 'sepbin':300, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='Hydra-A':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('Hydra-A_co21_22km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'Hydra-A_co21_22km_mom1.fits', 'ra':'09:18:05.68', 
                           'dec':'−12:05:43.8', 'res':1.053, 'nx':301, 'ny':301,
                           'ylab':np.array([4,5,6]), 'xlab':np.array([2,3,4,5,6,7]), 
                           'msk_sz':1000000, 'rkpc':1, 'xpl':2, 'xpu':7.5, 'ypl':3, 'ypu':7, 
                           'cbfrac':0.034, 'one_third':123, 'half':62, 'vsfyl':5, 'vsfyu':800, 
                           'histulim':20, 'sepbin':500, 'n_bins':200, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='2A0335':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('2A0335_20km_c010_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'2A0335_20km_c010_mom1.fits', 'ra':'03:38:40.54612', 
                           'dec':'+09:58:12.1373', 'res':0.700, 'nx':291, 'ny':291,
                           'ylab':np.array([10,20,30]), 'xlab':np.array([0,10,20,30,40]), 
                           'msk_sz':50000, 'rkpc':5, 'xpl':0, 'xpu':40, 'ypl':0, 'ypu':35, 
                           'cbfrac':0.040, 'one_third':62, 'half':35, 'vsfyl':10, 'vsfyu':300, 
                           'histulim':80, 'sepbin':300, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='A1795':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('a1795_10km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'a1795_10km_mom1.fits', 'ra':'13:48:52.495', 
                           'dec':'+26:35:34.32', 'res':1.220, 'nx':201, 'ny':201,
                           'ylab':np.array([0,10,20]), 'xlab':np.array([0,10,20]), 
                           'msk_sz':1000000, 'rkpc':8, 'xpl':0, 'xpu':25, 'ypl':0, 'ypu':25, 
                           'cbfrac':0.040, 'one_third':83, 'half':63, 'vsfyl':6, 'vsfyu':300, 
                           'histulim':80, 'sepbin':200, 'n_bins':100, 'cd':wcs.wcs.cd[1][1]*3600}
            # if pltmap=='flux':
        if gname=='A3581':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('a3581_20km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'a3581_20km_mom1.fits', 'ra':'14:07:29.76313', 
                           'dec':'-27:01:04.5689', 'res':0.435, 'nx':362, 'ny':361,
                           'ylab':np.array([0,5,10]), 'xlab':np.array([0,5,10]), 
                           'msk_sz':1000000, 'rkpc':3, 'xpl':0, 'xpu':10, 'ypl':0, 'ypu':10, 
                           'cbfrac':0.040, 'one_third':74, 'half':50, 'vsfyl':6, 'vsfyu':300, 
                           'histulim':80, 'sepbin':200, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':                
        if gname=='PKS0745':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('PKS0745-co10_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'PKS0745-co10_mom1.fits', 'ra':'07:47:31.321', 
                           'dec':'-19:17:39.97', 'res':1.890, 'nx':512, 'ny':512,
                           'ylab':np.array([100,105,110]), 'xlab':np.array([90,100,110,120]), 
                           'msk_sz':1000000, 'rkpc':8, 'xpl':90, 'xpu':120, 'ypl':95, 'ypu':115, 
                           'cbfrac':0.038, 'one_third':30, 'half':20, 'vsfyl':5, 'vsfyu':250, 
                           'histulim':80, 'sepbin':200, 'n_bins':100, 'cd':wcs.wcs.cd[1][1]*3600}
            # if pltmap=='flux':
        if gname=='R0821':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('RXJ0821_co10_20km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'RXJ0821_co10_20km_mom1.fits', 'ra':'08:21:02.25', 
                           'dec':'+07:51:47.4', 'res':2.000, 'nx':801, 'ny':801,
                           'ylab':np.array([10,20,30,40]), 'xlab':np.array([0,10,20,30,40]), 
                           'msk_sz':50000, 'rkpc':5, 'xpl':0, 'xpu':40, 'ypl':0, 'ypu':40, 
                           'cbfrac':0.033, 'one_third':15, 'half':8, 'vsfyl':2, 'vsfyu':200, 
                           'histulim':80, 'sepbin':500, 'n_bins':100, 'cd':wcs.wcs.cd[1][1]*3600}
            # if pltmap=='flux':
        if gname=='R1539':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('rxj1539_20km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'rxj1539_20km_mom1.fits', 'ra':'15:39:34.47', 
                           'dec':'-83:35:22.1', 'res':1.437, 'nx':116, 'ny':116,
                           'ylab':np.array([10,20,30,40,50]), 'xlab':np.array([0,20,40,60]), 
                           'msk_sz':80000, 'rkpc':8, 'xpl':0, 'xpu':50, 'ypl':0, 'ypu':60, 
                           'cbfrac':0.033, 'one_third':35, 'half':28, 'vsfyl':10, 'vsfyu':200, 
                           'histulim':100, 'sepbin':100, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='S1101':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('S1101_30km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'S1101_30km_mom1.fits', 'ra':'23:13:58.65', 
                           'dec':'-42:43:39.6', 'res':1.094, 'nx':362, 'ny':364,
                           'ylab':np.array([40,50]), 'xlab':np.array([30,40,50,60]), 
                           'msk_sz':1000000, 'rkpc':5, 'xpl':30, 'xpu':60, 'ypl':35, 'ypu':55, 
                           'cbfrac':0.040, 'one_third':55, 'half':35, 'vsfyl':7, 'vsfyu':200, 
                           'histulim':80, 'sepbin':100, 'n_bins':50, 'cd':wcs.wcs.cd[1][1]*3600}
            # if pltmap=='flux':
        if gname=='A1835': # nothing concrete on the resolution yet, no CDELT or CD
            if pltmap=='velocity':
                filename = get_pkg_data_filename('a1835-co10_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'a1835-co10_mom1.fits', 'ra':'14:01:02.08791', 
                           'dec':'+02:52:42.5689', 'res':3.966, 'nx':256, 'ny':256,
                           'ylab':np.array([160,170,180,190,200]), 'xlab':np.array([140,150,160,170,180]), 
                           'msk_sz':50000, 'rkpc':5, 'xpl':140, 'xpu':180, 'ypl':155, 'ypu':200, 
                           'cbfrac':0.040, 'one_third':45, 'half':25, 'vsfyl':10, 'vsfyu':200, 
                           'histulim':500, 'sepbin':200, 'n_bins':50, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='A1664':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('a1664-co10_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'a1664-co10_mom1.fits', 'ra':'09:18:05.68', 
                           'dec':'−12:05:43.8', 'res':2.302, 'nx':256, 'ny':256,
                           'ylab':np.array([20,40,60,80,100,120,140]), 'xlab':np.array([0,20,40,60,80,100,120,140]), 
                           'msk_sz':1000000, 'rkpc':1, 'xpl':0, 'xpu':150, 'ypl':0, 'ypu':150, 
                           'cbfrac':0.035, 'one_third':45, 'half':30, 'vsfyl':10, 'vsfyu':200, 
                           'histulim':250, 'sepbin':200, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='A262':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('a262_20km_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'a262_20km_mom1.fits', 'ra':'01:52:46.47', 
                           'dec':'+36:09:06.3', 'res':0.330, 'nx':231, 'ny':231,
                           'ylab':np.array([1,2,3,4,5,6,7]), 'xlab':np.array([0,1,2,3,4,5,6,7]), 
                           'msk_sz':1000000, 'rkpc':1, 'xpl':0, 'xpu':7, 'ypl':1, 'ypu':7.5, 
                           'cbfrac':0.035, 'one_third':102, 'half':52, 'vsfyl':7, 'vsfyu':500, 
                           'histulim':20, 'sepbin':500, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
        if gname=='Phoenix-A':
            if pltmap=='velocity':
                filename = get_pkg_data_filename('phoenix_30_mom1.fits')
                hdu = fits.open(filename)[0]
                wcs = WCS(hdu.header)
                sys_par = {'vmapfn':'phoenix_30_mom1.fits', 'ra':'23:44:43.94576', 
                           'dec':'-42:43:12.2950', 'res':6.750, 'nx':362, 'ny':364,
                           'ylab':np.array([5,10,15]), 'xlab':np.array([0,5,10,15]), 
                           'msk_sz':1000000, 'rkpc':5, 'xpl':0, 'xpu':15, 'ypl':0, 'ypu':15, 
                           'cbfrac':0.040, 'one_third':80, 'half':55, 'vsfyl':10, 'vsfyu':300, 
                           'histulim':20, 'sepbin':200, 'n_bins':100, 'cd':wcs.wcs.cdelt[1]*3600}
            # if pltmap=='flux':
    return sys_par
