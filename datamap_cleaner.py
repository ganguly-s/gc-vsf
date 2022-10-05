import numpy as np
import os
import matplotlib.pyplot as plt
import datamap_mask as mask

def make_data(alldata):
    data, error_data = alldata['velo_data'], alldata['error_data']
    good_v = alldata['good_v']
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
    rB_half = np.ndarray.flatten(np.triu(rB_matrix, k=0))
    rA=rA_half[good_dist]
    rB=rB_half[good_dist]
    inner_f=(rA<=r)&(rB<=r)
    outer_f=(rA>r)&(rB>r)
    np_dist_in=np_dist[inner_f]
    np_v_diff_in=np_v_diff[inner_f]
    np_error_in_2=np_error_2[inner_f]
    np_dist_out=np_dist[outer_f]
    np_v_diff_out=np.abs(np_v_diff[outer_f])
    np_error_out_2=np_error_2[outer_f]
    return np_dist, np_v_diff, np_error_2, np_dist_in, np_v_diff_in, np_error_in_2, np_dist_out, np_v_diff_out, np_error_out_2

def calc_vsf(alldata):
    n_bins = 200
    gname, telescope = alldata['gname'], alldata['telescope']
    impath = telescope+'/vsfplots/'+gname+'/'
    datpath = telescope+'/vsfdat/'+gname+'/'
    if not os.path.exists(telescope+'/'):
        os.mkdir(telescope+'/')
    if not os.path.exists(telescope+'/vsfplots/'):
        os.mkdir(telescope+'/vsfplots/')
    if not os.path.exists(telescope+'/vsfplots/'+gname+'/'):
        os.mkdir(telescope+'/vsfplots/'+gname+'/')
    fnvsfall = impath+gname+'_allvsf.png'
    if not os.path.exists(telescope+'/'):
        os.mkdir(telescope+'/')
    if not os.path.exists(telescope+'/vsfdat/'):
        os.mkdir(telescope+'/vsfdat/')
    if not os.path.exists(telescope+'/vsfdat/'+gname+'/'):
        os.mkdir(telescope+'/vsfdat/'+gname+'/')
    fnvsfallnpz = datpath+gname+'_allvsf'
    rkpc, res = alldata['rkpc'], alldata['res']
    r = rkpc/res
    if 'vsfbin' in alldata:
        n_bins = alldata['vsfbin']
    one_third, half = alldata['one_third'], alldata['half']
    vsfyl, vsfyu, vsfbin = alldata['vsfyl'], alldata['vsfyu'], alldata['vsfbin']
    n_bins_in = 50
    n_bins_out = 50
    (np_dist, np_v_diff, np_error_2, np_dist_in, np_v_diff_in, np_error_in_2, np_dist_out, np_v_diff_out, np_error_out_2) = make_data(alldata)    
    d_max = np.max(np_dist)
    d_max_in = np.max(np_dist_in)
    d_max_out = np.max(np_dist_out)
    # need to be careful with creating the list of l. The first few values need to be set by hand because of the discretization of the map.
    dist_floor = np.floor(np_dist*100)
    unique = np.unique(dist_floor)/100.
    dist_array = np.append(unique[0:15], np.logspace(np.log10(unique[15]), np.log10(d_max), n_bins-15))
    dist_array_in = np.append(unique[0:15], np.logspace(np.log10(unique[15]), np.log10(d_max_in), n_bins_in-15))
    dist_array_out = np.append(unique[0:15], np.logspace(np.log10(unique[15]), np.log10(d_max_out), n_bins_out-15))
    
    dist_array_kpc = dist_array*res
    dist_array_kpc_in = dist_array_in*res
    dist_array_kpc_out = dist_array_out*res
    y_expect = dist_array**(1.0/3)*one_third
    y_expect2 = dist_array**(1.0/2)*half
    y_expect3 = dist_array**(2./3)*14
    y_expect4 = dist_array*4
    
    v_diff_mean = np.zeros(n_bins)
    v_diff_mean2 = np.zeros(n_bins)
    error_mean = np.zeros(n_bins)
    v_diff_mean_in=np.zeros(n_bins_in)
    v_diff_mean_out=np.zeros(n_bins_out)
    error_mean_in=np.zeros(n_bins_in)
    error_mean_out=np.zeros(n_bins_out)
    for i in range(0, n_bins-1):
        this_bin = (np_dist >= dist_array[i]) & (np_dist < dist_array[i+1])
        v_diff_mean[i] = np.mean(np.abs(np_v_diff[this_bin]))
        v_diff_mean2[i] = np.mean((np_v_diff[this_bin])**2)
        error_mean[i] = np.sqrt(np.sum(np_error_2[this_bin]))/np_error_2[this_bin].size#np.sqrt(np.mean(np_error_2[this_bin]))
    for i in range(0, n_bins_in-1):
        this_bin=(np_dist_in>=dist_array_in[i])&(np_dist_in<dist_array_in[i+1])
        v_diff_mean_in[i]=np.mean(np.abs(np_v_diff_in[this_bin]))
        error_mean_in[i] = np.sqrt(np.sum(np_error_in_2[this_bin]))/np_error_in_2[this_bin].size#np.sqrt(np.mean(np_error_in_2[this_bin]))
    for i in range(0, n_bins_out-1):
        this_bin=(np_dist_out>=dist_array_out[i])&(np_dist_out<dist_array_out[i+1])
        v_diff_mean_out[i]=np.mean(np.abs(np_v_diff_out[this_bin]))
        error_mean_out[i] = np.sqrt(np.sum(np_error_out_2[this_bin]))/np_error_out_2[this_bin].size#np.sqrt(np.mean(np_error_out_2[this_bin]))
    error_mean_smooth = mask.run_helper(dist_array_kpc, error_mean)
    v_diff_mean_smooth = mask.run_helper(dist_array_kpc, v_diff_mean)
    lower_v = v_diff_mean_smooth-error_mean_smooth
    upper_v = v_diff_mean_smooth+error_mean_smooth
    error_mean_smooth_in = mask.run_helper(dist_array_kpc_in, error_mean_in)
    v_diff_mean_smooth_in = mask.run_helper(dist_array_kpc_in, v_diff_mean_in)
    lower_v_in = v_diff_mean_smooth_in-error_mean_smooth_in
    upper_v_in = v_diff_mean_smooth_in+error_mean_smooth_in
    error_mean_smooth_out = mask.run_helper(dist_array_kpc_out, error_mean_out)
    v_diff_mean_smooth_out = mask.run_helper(dist_array_kpc_out, v_diff_mean_out)
    lower_v_out = v_diff_mean_smooth_out-error_mean_smooth_out
    upper_v_out = v_diff_mean_smooth_out+error_mean_smooth_out
    np.savez(fnvsfallnpz, rkpc = rkpc, dist_array = dist_array, dist_array_kpc=dist_array_kpc, dist_array_kpc_in=dist_array_kpc_in, dist_array_kpc_out=dist_array_kpc_out,
                 v_diff_mean_smooth=v_diff_mean_smooth, v_diff_mean_smooth_in=v_diff_mean_smooth_in, upper_v = upper_v, upper_v_in = upper_v_in, upper_v_out = upper_v_out,
                 v_diff_mean_smooth_out=v_diff_mean_smooth_out, lower_v = lower_v, lower_v_in = lower_v_in, lower_v_out = lower_v_out)
    return