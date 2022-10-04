import numpy as np

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