#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:12:52 2022

@author: sg
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import datamap_mask as mask
import inputs as inp

newparams = {'font.family':'serif', 'axes.labelsize': 24, 'axes.linewidth': 1, 
              'lines.linewidth': 2, 'figure.figsize': (8, 6),
              'figure.subplot.wspace': 0.4,
              'ytick.labelsize': 24, 'xtick.labelsize': 24,
              'ytick.major.pad': 5, 'xtick.major.pad': 5,
              'legend.fontsize': 24, 'legend.frameon': True,
              'legend.handlelength': 2, 'figure.facecolor':'white',
              'savefig.dpi': 300}
plt.rcParams.update(newparams)

gname = 'S1101'
telescope = 'MUSE'
# telescope = 'ALMA'
pltmap = 'velocity'
# pltmap = 'flux'
max_error = 5
sysparam = inp.get_sys_params(gname,telescope,pltmap)
one_third, half = sysparam['one_third'], sysparam['half']
vsfyl, vsfyu, vsfbin = sysparam['vsfyl'], sysparam['vsfyu'], sysparam['vsfbin']

gname = ['2A0335','Centaurus','A1795','A3581','Hydra-A','PKS0745','R0821','R1539','S1101']
colors = plt.cm.jet(np.linspace(0,1,len(gname)))
#impath = telescope+'/vsfplots/'+gname+'/'
impath = '../final_images/'
#datpath = telescope+'/vsfdat/'+gname+'/'
if not os.path.exists(impath):
    os.mkdir(impath)
#if not os.path.exists(telescope+'/vsfplots/'):
#    os.mkdir(telescope+'/vsfplots/')
#if not os.path.exists(telescope+'/vsfplots/'+gname+'/'):
#    os.mkdir(telescope+'/vsfplots/'+gname+'/')
#fnvsfall = impath+gname+'_allvsf.png'
fnvsfall = impath+'fig1.pdf'
#fnvsfallnpz = datpath+gname+'_allvsf.npz'
#fnvsfallnpzunsmooth = datpath+gname+'_allvsf_unsmooth.npz'

plt.figure(figsize=(13, 12))
y_expect = np.array([0.01,1,150])**(1.0/3)*150#one_third
y_expect2 = np.array([0.01,1,150])**(1.0/2)*150#half

plt.loglog([0.01,1,150], y_expect, "k--", label="1/3")
plt.loglog([0.01,1,150], y_expect2, "m-.", label="1/2")

for i,gn in enumerate(gname):
    datpath = telescope+'/vsfdat/'+gn+'/'
    fnvsfallnpz = datpath+gn+'_allvsf.npz'
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

    plt.loglog(dist_array_kpc, v_diff_mean_smooth, marker="o",color=colors[i],linestyle="None", markersize=4, label=gn)
    plt.fill_between(dist_array_kpc[:-1], lower_v[:-1], upper_v[:-1],color=colors[i], alpha=0.5)

# plt.loglog(dist_array_kpc_in, v_diff_mean_smooth_in, marker="s",linestyle="None", markersize=4, color="C3",label='inner (r<'+str(rkpc)+' kpc)')
# plt.fill_between(dist_array_kpc_in[:-1], lower_v_in[:-1], upper_v_in[:-1], color="C3", alpha=0.2)
# plt.loglog(dist_array_kpc_out, v_diff_mean_smooth_out, marker="D",linestyle="None", markersize=4, color="C2",label='outer (r>'+str(rkpc)+' kpc)')
# plt.fill_between(dist_array_kpc_out[:-1], lower_v_out[:-1], upper_v_out[:-1], color="C2", alpha=0.2)

#vsfdata= np.load(fnvsfallnpzunsmooth)
#rkpc = vsfdata['rkpc']
#dist_array = vsfdata['dist_array']
#dist_array_kpc = vsfdata['dist_array_kpc']
#dist_array_kpc_in = vsfdata['dist_array_kpc_in']
#dist_array_kpc_out = vsfdata['dist_array_kpc_out']
#v_diff_mean_unsmooth = vsfdata['v_diff_mean_smooth']
#v_diff_mean_smooth_in = vsfdata['v_diff_mean_smooth_in']
#v_diff_mean_smooth_out = vsfdata['v_diff_mean_smooth_out']
#lower_v = vsfdata['lower_v']
#lower_v_in = vsfdata['lower_v_in']
#lower_v_out = vsfdata['lower_v_out']
#upper_v = vsfdata['upper_v']
#upper_v_in = vsfdata['upper_v_in']
#upper_v_out = vsfdata['upper_v_out']

#plt.loglog(dist_array_kpc, v_diff_mean_unsmooth, marker="o",linestyle="None", markersize=4, color="C1",label=gname+' All (unsmooth)')
#plt.fill_between(dist_array_kpc[:-1], lower_v[:-1], upper_v[:-1], color="C1", alpha=0.2)

plt.xlabel("separation (kpc)")
plt.ylabel(r"$|\delta v|\, \rm (km/s)$")
plt.legend(loc="center left", prop={'size': 22}, bbox_to_anchor=(1, 0.5))
#plt.ylim(vsfyl, vsfyu)
plt.ylim(5,350)
plt.xlim(0.03,150)
plt.grid()
plt.tick_params(which='both',direction='in')
#plt.title(gname, size=24)
plt.savefig(fnvsfall,bbox_inches="tight")
plt.show()
