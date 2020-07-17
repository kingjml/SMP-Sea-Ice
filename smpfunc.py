import numpy as np
import pandas as pd
import math
from scipy.signal import resample
from scipy.ndimage import gaussian_filter, interpolation

def scale_profile(scaling, depth_array, value_array, l_resample = 40, h_resample = 1):
    """Scale snow property profiles to new distance representation
    Argrumnets: 
        scaling: An array of % values to scale each layer by
        depth_array: The origional distance array
        value_array: The origoinal snow property array
        l_resample: Height of layers in mm
        h_resample: Resampled resolution output in mm
    Output:
        result_dist: Resampled and scaled distance
        result_value: REsampled and scaled snow property
    """
    result_value = np.array([])
    result_dist = np.array([])
    
    delta_h = np.diff(depth_array)[0] #This assumes its equal!
    if h_resample is None:
        h_resample = delta_h
    l_last = np.array(-h_resample)

    layer_thickness = l_resample + (scaling * l_resample)
    new_total_thickness = sum(layer_thickness)
    total_stretch = abs((new_total_thickness - depth_array.max())/depth_array.max())
    depth_stretch= np.arange(0,new_total_thickness, h_resample)

    for l_idx in np.arange(len(layer_thickness)):
        ol_start = l_resample*l_idx
        old_end = ol_start + l_resample
        ol_start_idx = (np.abs(depth_array - ol_start)).argmin() # Do this incase mod(delta_h, l_resample)!=0 
        ol_end_idx = (np.abs(depth_array - old_end)).argmin()+1

        orig_dist = depth_array[ol_start_idx:ol_end_idx]
        orig_value = value_array[ol_start_idx:ol_end_idx]
        
        l_thickness = np.floor(layer_thickness[l_idx])
        l_end = np.floor(layer_thickness[0:l_idx+1].sum())
        l_end_idx = (np.abs(depth_stretch - l_end)).argmin()
        l_start = l_last + h_resample
        l_start_idx = (np.abs(depth_stretch - l_start)).argmin()
        l_last = l_end
        
        if len(orig_value) > 2 :
            scaled_dist = depth_stretch[l_start_idx:l_end_idx+1]
            scaled_value, x1 = resample(orig_value,len(scaled_dist),orig_dist) 
            result_dist = np.append(result_dist,scaled_dist)
            result_value = np.append(result_value, scaled_value) 
    
    return result_dist, result_value

def score_scaling(scaling, depth_array, value_array, obs_array, l_resample = 40, h_resample = 1):
    result_dist, result_rho = scale_profile(scaling, depth_array, value_array, l_resample, h_resample)
    smp_rho_interp = np.interp(depth_array,result_dist,result_rho)
    return (1/result_dist.max()*(np.power(np.log10(smp_rho_interp) - np.log10(obs_array),2)).sum())

def random_stretch(layer_num, max_change = 0.1, max_change_layer = 0.7):
    """Generate random numbers to manipulate either stretch (+) or erode
    Argrumnets: 
        layer_num: Number of random values to generate, ie. number of snowpack layers
        max_change: The maximum precentage that the total snowpack can be modified
        max_change_layer: The maximum precentage that an individual layer can be modified
    Output:
        A list of modifiers to erode or dilate snowpack layers based on the inputs
    """
    while True:
        layer_stretch = np.array([])
        layer_order = np.random.choice(layer_num, layer_num, replace=False)
        for layer in np.arange(0,len(layer_order)):
            layer_stretch =  np.append(layer_stretch, np.random.uniform(-max_change_layer,max_change_layer))
        if (np.abs(layer_stretch.sum()) <= max_change):
            break
    layer_order_sort = np.argsort(layer_order)
    return layer_stretch[layer_order_sort]

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return np.isnan(y), lambda z: z.to_numpy().nonzero()[0]

def preprocess(profile, smoothing = 0.5, noise_threshold = 0.01, verbose = False):
    # Check for values below 0.01N, remove them from the analysis
    if verbose:
        print('Values below 0.01N removed: {}'.format((profile.samples.force < noise_threshold).sum()))
    num = profile.samples.force._get_numeric_data()
    profile.samples['force'] = gaussian_filter(profile.samples['force'] , sigma=smoothing)
    num[num < noise_threshold] = np.nan
    nans, x= nan_helper(profile.samples.force)
    profile.samples.force[nans]= np.interp(x(nans), x(~nans), profile.samples.force[~nans])
    return profile

def rolling_window(a, window, fun=None, pad=False):
    window = (np.ceil(window) // 2 * 2 + 1).astype(int) #round up to next odd number
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    rWindow = np.lib.stride_tricks.as_strided(a, shape=shape, strides=a.strides + (a.strides[-1],))
    if fun:
        rWindow = fun(rWindow, -1)
    if pad: #This will be slow if no function is applied!
      padSize = int(np.absolute(rWindow.shape[0]-a.shape[0])/2)
      rWindow = np.lib.pad(rWindow, (padSize,padSize), 'constant', constant_values=np.nan)
    return rWindow

def extract_samples(a_height, a_value, b_height, window_size):
    return pd.DataFrame({'count_samp': [np.count_nonzero(a_value[np.abs(a_height-height)<=window_size]) for height in b_height],
                         'mean_samp': [np.mean(a_value[np.abs(a_height-height)<=window_size]) for height in b_height], 
                         'median_samp': [np.median(a_value[np.abs(a_height-height)<=window_size]) for height in b_height],
                         'stdev_samp': [np.std(a_value[np.abs(a_height-height)<=window_size]) for height in b_height]})

def calc_skill(result, window_size, drop_na = True):
    if drop_na:
        result.dropna(inplace=True)
    result = result[result['count_samp']>=window_size] # Remove comparisons outside the profile
    result = result[~result['TYPE'].isin(['N', 'I'])] # Remote new snow and ice because we don't have enough samples
    
    r = np.corrcoef(result['mean_samp'],result['RHO'])[1][0]
    rmse = np.sqrt(np.mean((result['mean_samp']-result['RHO'])**2))
    rmse_corr = rmse - np.abs(np.mean(result['mean_samp']-result['RHO']))
    mae = np.sum(np.abs(result['mean_samp']-result['RHO']))/np.ma.count(result['RHO'])
    
    return r, rmse, rmse_corr, mae