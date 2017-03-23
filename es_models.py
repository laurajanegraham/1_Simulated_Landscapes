# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:10:21 2017

@author: lg1u16
"""
import sys
sys.path.insert(0, 'SCALEFORES/1_Simulated_Landscapes')

import nlmfunctions as nlm
import pandas as pd

# this section sets up all of the information about the NLM created
nRow = 51 # Number of rows
nCol = 51 # Number of columns
p = range(1, 10, 1) # range of values of p (proportion of NLC)
h = range(0, 11, 1) # range of values of h (spatial autocorrelation)
r = range(100) # number of replicates per landscape type
w = [3,5,9,15] # window sizes for moving window analysis
sigma = [0.02, 0.03, 0.06, 0.12, 0.5, 0.75, 1]
ex_keys = [{'lc':1, 'a': 0.0000001}, {'lc':1, 'a':1}, {'lc':1, 'a':5}, {'lc':1, 'a':50}, {'lc':1, 'a':-1}, {'lc':1, 'a':-5}, {'lc':1, 'a':-50}]

# create the parameters to create surfaces for based on how they were set up above
param_set = pd.DataFrame(nlm.expandgrid(p, h, r, w, ex_keys))
param_set = param_set.rename(index=str, columns = {'Var1':'p', 'Var2':'h', 'Var3':'r', 'Var4':'w', 'Var5':'ex_keys'})
param_set['p'] = param_set['p']/10
param_set['h'] = param_set['h']/10

## run the model on the full set of parameters indiscriminate of focal patch land-cover type (recreation?)
#es_mod_rec = param_set.apply(nlm.simple_esmod, args=(nRow, nCol), axis=1)
#es_mod_rec = es_mod_rec.rename(index=str, columns = {'p': 'p_val', 'h': 'h_val', 'r': 'rep', 'w': 'window_size', 'es_mean': 'es_mean', 'es_total': 'es_total', 'es_var': 'es_var'})
#es_mod_rec.to_csv('C:/Users/lg1u16/SCALEFORES/1_Simulated_Landscapes/results/es_mod_rec.csv', index=False)
#
## run the model where focal patch has to be natural (==1) (e.g. pollination)
#es_mod_poll_agri = param_set.apply(nlm.two_step_esmod, args=(nRow, nCol), axis=1)
#es_mod_poll_agri = es_mod_poll_agri.rename(index=str, columns = {'p': 'p_val', 'h': 'h_val', 'r': 'rep', 'w': 'window_size', 'es1_mean': 'es1_mean', 'es1_total': 'es1_total', 'es1_var': 'es1_var', 'es2_mean': 'es2_mean', 'es2_total': 'es2_total', 'es2_var': 'es2_var'})
#es_mod_poll_agri.to_csv('C:/Users/lg1u16/SCALEFORES/1_Simulated_Landscapes/results/es_mod_poll_agri.csv', index=False)

# experiment with changing the function in the simple model
es_mod_rec2 = param_set.apply(nlm.simple_esmod, args=(nlm.exp_func, nRow, nCol), axis=1)
es_mod_rec2 = es_mod_rec2.rename(index=str, columns = {'p': 'p_val', 'h': 'h_val', 'r': 'rep', 'w': 'window_size', 'es_mean': 'es_mean', 'es_total': 'es_total', 'es_var': 'es_var'})
es_mod_rec2.to_csv('C:/Users/lg1u16/SCALEFORES/1_Simulated_Landscapes/results/es_mod_rec_exp.csv', index=False)

# experiment with applying a distance decay function 
dd_param_set = pd.DataFrame(nlm.expandgrid(p, h, r, sigma))
dd_param_set = dd_param_set.rename(index=str, columns = {'Var1':'p', 'Var2':'h', 'Var3':'r', 'Var4':'sigma'})
dd_param_set['p'] = dd_param_set['p']/10
dd_param_set['h'] = dd_param_set['h']/10

es_mod_dd = dd_param_set.apply(nlm.dd_simple_esmod, args=(nRow, nCol), axis=1)
es_mod_dd = es_mod_dd.rename(index=str, columns =  {'p': 'p_val', 'h': 'h_val', 'r': 'rep', 'sigma': 'sigma', 'es1_mean': 'es1_mean', 'es1_total': 'es1_total', 'es1_var': 'es1_var', 'es2_mean': 'es2_mean', 'es2_total': 'es2_total', 'es2_var': 'es2_var'})
es_mod_dd.to_csv('C:/Users/lg1u16/SCALEFORES/1_Simulated_Landscapes/results/es_mod_dd.csv', index=False)

# get some example landscapes
param_set_sml = param_set.groupby(['p', 'h']).size().reset_index()
es_mod_example = pd.DataFrame()
p_range = [0.1, 0.5, 0.9]
h_range = [0, 0.5, 0.9]
for p in p_range:
    for h in h_range:
        es_mod_example = es_mod_example.append(mpd_prop_df(25, 25, p, h))

es_mod_example.to_csv('SCALEFORES/1_Simulated_Landscapes/es_mod_example.csv')