# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 11:17:32 2017

@author: lg1u16
"""

from nlmpy import nlmpy # for some reason my installation requires me to load it like this...
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from scipy.ndimage.filters import generic_filter

# this section sets up all of the information about the NLM created
nRow = 50 # Number of rows
nCol = 50 # Number of columns
p = range(1, 10, 1) # range of values of p (proportion of NLC)
h = range(0, 11, 1) # range of values of h (spatial autocorrelation)
r = range(100) # number of replicates per landscape type
w = [1,3,5,9,15] # window sizes for moving window analysis

# SET UP FUNCTIONS #############################################################
# function to replicate expand.grid from R
def expandgrid(*itrs):
   product = list(itertools.product(*itrs))                                      
   return {'Var{}'.format(i+1):[x[i] for x in product] for i in range(len(itrs))}

# function to get binary classification from MPDM
def mpd_prop(nRow, nCol, h, p):
    mpd_out = nlmpy.mpd(nRow, nCol, h)
    prop_out = nlmpy.classifyArray(mpd_out, [1 - p, p])
    return prop_out;
    
# function to get binary classification from MPDM - put into df form
def mpd_prop_df(p, h):
    mpd_out = nlmpy.mpd(250, 250, h)
    prop_out = nlmpy.classifyArray(mpd_out, [1 - p, p])
    out = pd.DataFrame(prop_out)
    out['y'] = range(0, nRow, 1)
    out = pd.melt(out, id_vars=['y'], var_name='x', value_name='value')
    out['h'] = h
    out['p'] = p
    return out;

# function to bring together creating the landscape and moving window analysis
def simple_esmod(params, nRow, nCol):
    p = params['p']
    h = params['h']
    w = int(params['w'])
    r = params['r']
    out = mpd_prop(nRow, nCol, h, p)

    wdw = generic_filter(out, np.mean, w, mode='wrap')
    # this is currently set to take the relationship as 1:1 (i.e. 10% natural cover = 10% ecosystem service)
    # will need to add in the relationship to create ES surface at a later date
    es_mean = np.mean(wdw)
    es_total = np.sum(wdw)
    es_var = np.var(wdw) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1
    return pd.Series({'p_val': p, 'h_val': h, 'rep': r, 'window_size': w, 'es_mean': es_mean, 'es_total': es_total, 'es_var': es_var});

# function to bring together creating the landscape and moving window analysis in a two step approach (e.g. pollination as input to agri)
def two_step_esmod(params, nRow, nCol):
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
 
# create the parameters to create surfaces for based on how they were set up above
param_set = pd.DataFrame(expandgrid(p, h, r, w))
param_set = param_set.rename(index=str, columns = {'Var1':'p', 'Var2':'h', 'Var3':'r', 'Var4':'w'})
param_set['p'] = param_set['p']/10
param_set['h'] = param_set['h']/10

# run the model on the full set of parameters indiscriminate of focal patch land-cover type (recreation?)
es_mod_rec = param_set.apply(simple_esmod, args=(nRow, nCol), axis=1)
es_mod_rec = es_mod_rec.rename(index=str, columns = {'p': 'p_val', 'h': 'h_val', 'r': 'rep', 'w': 'window_size', 'es_mean': 'es_mean', 'es_total': 'es_total', 'es_var': 'es_var'})
es_mod_rec.to_csv('C:/Users/lg1u16/SCALEFORES/1_Simulated_Landscapes/results/es_mod_rec.csv', index=False)

# run the model where focal patch has to be natural (==1) (e.g. pollination)
es_mod_poll_agri = param_set.apply(two_step_esmod, args=(nRow, nCol), axis=1)
es_mod_poll_agri = es_mod_poll_agri.rename(index=str, columns = {'p': 'p_val', 'h': 'h_val', 'r': 'rep', 'w': 'window_size', 'es1_mean': 'es1_mean', 'es1_total': 'es1_total', 'es1_var': 'es1_var', 'es2_mean': 'es2_mean', 'es2_total': 'es2_total', 'es2_var': 'es2_var'})
es_mod_poll_agri.to_csv('C:/Users/lg1u16/SCALEFORES/1_Simulated_Landscapes/results/es_mod_poll_agri.csv', index=False)

# get some example landscapes
#param_set_sml = param_set.groupby(['p', 'h']).size().reset_index()
#es_mod_example = pd.DataFrame()
#for p in range(1, 10):
#    p = p/10
#    for h in range(0, 11):
#        h = h/10
#        es_mod_example = es_mod_example.append(mpd_prop_df(p, h))
#
#es_mod_example.to_csv('SCALEFORES/1_Simulated_Landscapes/es_mod_example.csv')