# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 11:17:32 2017

@author: lg1u16
"""

from nlmpy import nlmpy # for some reason my installation requires me to load it like this...
import numpy as np
import pandas as pd
import itertools
from scipy.ndimage.filters import generic_filter
import math


# SET UP FUNCTIONS #############################################################
# function to replicate expand.grid from R
def expandgrid(*itrs):
    """Copy the R expandgrid function. Takes 2+ variables and outputs the unique 
    combinations of these
    
    Parameters
    ----------
    *itrs: vectors
        The variables to be included
    
    Returns
    -------
    expanded_grid: array
        Each row is a unique combination of the supplied vectors
    """
    product = list(itertools.product(*itrs))    
    expandedgrid = {'Var{}'.format(i+1):[x[i] for x in product] for i in range(len(itrs))}                                  
    return expandedgrid

# function to get binary classification from MPDM
def mpd_prop(nRow, nCol, h, p):
    """Create a neutral landscape which has a specific spatial autocorrelation 
    and proportion
    
    Parameters
    ----------
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    h: float
        Level of spatial autocorrelation - 0 is random, 1 is totally correlated
    p: float
        Proportion of landcover 1 in the landscape
        
    Returns
    -------
    prop_out: np.array
        NLM with required features
    """
    mpd_out = nlmpy.mpd(nRow, nCol, h)
    prop_out = nlmpy.classifyArray(mpd_out, [1 - p, p])
    return prop_out;
    
# function to get binary classification from MPDM - put into df form
def mpd_prop_df(nRow, nCol, p, h):
    """Create a neutral landscape which has a specific spatial autocorrelation 
    and proportion. This outputs as a dataframe for plotting in R.
    
    Parameters
    ----------
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    h: float
        Level of spatial autocorrelation - 0 is random, 1 is totally correlated
    p: float
        Proportion of landcover 1 in the landscape
        
    Returns
    -------
    out: pd.DataFrame
        NLM with required features
    """
    
    mpd_out = nlmpy.mpd(nRow, nCol, h)
    prop_out = nlmpy.classifyArray(mpd_out, [1 - p, p])
    out = pd.DataFrame(prop_out)
    out['y'] = range(0, nRow, 1)
    out = pd.melt(out, id_vars=['y'], var_name='x', value_name='value')
    out['h'] = h
    out['p'] = p
    return out;

# function to bring together creating the landscape and moving window analysis
def simple_esmod(params, uFunction, nRow, nCol):
    """Simplest ES model where the ES value is the value of the mean land cover 
    within a particular window. A function can be applied to this using the 
    uFunction argument.
    
    Parameters
    ----------
    params: pd.series
        Values of each of the parameters (p, h, w - window size, r - replicate, 
        and ex_keys - arguments to the user defined function)
    uFunction: function
        The function to be applied to calculate the ES value
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    """
    p = params['p']
    h = params['h']
    w = int(params['w'])
    r = params['r']
    ex_keys = params['ex_keys']
    
    out = mpd_prop(nRow, nCol, h, p)
    
    wdw = generic_filter(out, uFunction, w, mode='wrap', extra_keywords = ex_keys)
    
    # output values
    es_mean = np.mean(wdw)
    es_total = np.sum(wdw)
    es_var = np.var(wdw) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1
    return pd.Series({'p_val': p, 'h_val': h, 'rep': r, 'window_size': w, 'function': uFunction, 'ex_keys': ex_keys, 'es_mean': es_mean, 'es_total': es_total, 'es_var': es_var});

def dd_simple_esmod(params, nRow, nCol):
    """Simplest ES model where the ES value a distance weighting on the value of 
    the land cover within a particular window. ES1 is currently comparable to 
    the boring one and ES2 is comparable to ES1 in the two step model
    
    Parameters
    ----------
    params: pd.series
        Values of each of the parameters (p, h, sigma - the 'scale' parameter 
        for the distance decay function, r - replicate)
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    
    """
    p = params['p']
    h = params['h']
    r = params['r']
    sigma = params['sigma']

    # simulated landscape
    out = mpd_prop(nRow, nCol, h, p)
    
    # create the weighting surface which will be included as an input to the function
    x, y = np.meshgrid(np.arange(nRow), np.arange(nRow))
    x, y = (x/(nRow-1), y/(nRow-1))
    centre = (0.5, 0.5)
    #centre = (math.floor(nRow/2), math.floor(nRow/2))
    d = np.sqrt((x-centre[1])**2 + (y-centre[0])**2)
    wt = np.exp(-d**2/2*sigma**2)
    wt = wt/np.sum(wt)
    wt = wt.flatten()
      
    wdw1 = generic_filter(out, dd_func, nRow, mode='wrap', extra_keywords = {'lc':1, 'wt':wt})
    wdw2 = wdw1*out
    # output values
    es1_mean = np.mean(wdw1)
    es1_total = np.sum(wdw1)
    es1_var = np.var(wdw1) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1#
    
    es2_mean = np.mean(wdw2)
    es2_total = np.sum(wdw2)
    es2_var = np.var(wdw2)
    return pd.Series({'p_val': p, 'h_val': h, 'rep': r, 'sigma': sigma, 'es1_mean': es1_mean, 'es1_total': es1_total, 'es1_var': es1_var, 'es2_mean': es2_mean, 'es2_total': es2_total, 'es2_var': es2_var});

# function to bring together creating the landscape and moving window analysis in a two step approach (e.g. pollination as input to agri)
def two_step_esmod(params, nRow, nCol):
    """Two stage ES model where first an intermediary layer is calculated as the 
    same as the simple ES model but only if the land cover type is that required. 
    The ES value is then calculated as the mean of the intermediary layer within 
    the desired window.
    
    Parameters
    ----------
    params: pd.series
        Values of each of the parameters (p, h, w - window size, r - replicate)
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape    
    """
    p = params['p']
    h = params['h']
    w = int(params['w'])
    r = params['r']
    out = mpd_prop(nRow, nCol, h, p)
    
    # define the ES function
    # output the first layer (e.g. the one where focal patch must be = 1)
    wdw1 = generic_filter(out, np.mean, w, mode='wrap')
    # multiply wdw by out to set the zeros to zero
    wdw1 = wdw1 * out
    # this is currently set to take the relationship as 1:1 (i.e. 10% natural cover = 10% ecosystem service)
    # will need to add in the relationship to create ES surface at a later date
    es1_mean = np.mean(wdw1)
    es1_total = np.sum(wdw1)
    es1_var = np.var(wdw1) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1
    
    # output the second layer (e.g. the one which takes in the first as input and works on the opposite land cover = 0)
    wdw2 = generic_filter(wdw1, np.mean, w, mode='wrap')
    wdw2 = wdw2 * (1 - out)
    es2_mean = np.mean(wdw2)
    es2_total = np.sum(wdw2)
    es2_var = np.var(wdw2) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1
    
    return pd.Series({'p_val': p, 'h_val': h, 'rep': r, 'window_size': w, 'es1_mean': es1_mean, 'es1_total': es1_total, 'es1_var': es1_var, 'es2_mean': es2_mean, 'es2_total': es2_total, 'es2_var': es2_var});

# create functions for exploring the impact of the shape of the relationship between LC proportion and ES output
def lc_prop(values, lc):
    """Calculate the ES value based on a straight line relationship between the 
    mean LC value and the ES value.
    
    Parameters
    ----------
    values: array
        The array to process
    lc: int
        The land cover number to calculate the proportion of
    
    Returns
    -------
    es_value: float
        This is the value of the ES based on the defined rules {0, 1}  
    """
    out = int((values == lc).sum()) / values.size
    return out

def exp_func(values, lc, a):
    """Calculate the ES value based on an exponential relationship between the 
    mean LC value and the ES value. The returned value is scaled between 0 and 1.
    
    Parameters
    ----------
    values: array
        The array to process
    lc: int
        The land cover number to calculate the proportion of
    a: float
        The slope of the exponential function - the larger a is, the steeper the slope 
        a > 0 has exponential rise, a < 0 has exponential decline
        
    Returns
    -------
    es_value: float
        This is the value of the ES based on the defined rules  {0, 1}
    """
    lc_prop = int((values == lc).sum()) / values.size
    exp_value = np.exp(a*lc_prop)
    es_value = (exp_value - np.exp(a*0)) / (np.exp(a*1) - np.exp(a*0))
    return es_value

def dd_func(values, lc, wt):
    """Calculate the ES value where the desired LC in each cell is weighted by 
    it's distance to the central cell. 
    
    Parameters
    ----------
    values: array
        The array to process
    lc: int
        The land cover number to calculate the proportion of
    w: array
        A distance weighted surface
    
    Returns
    -------
    es_value: float
        This is the value of the ES based on the defined rules  {0, 1}
    """
    lc_one = (values == lc)*1
    es_value = (lc_one*wt).sum()
    return es_value
    
    
