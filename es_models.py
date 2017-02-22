# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:10:21 2017

@author: lg1u16
"""
import sys
sys.path.insert(0, 'SCALEFORES/1_Simulated_Landscapes')

import nlmfunctions as nlm

# this section sets up all of the information about the NLM created
nRow = 50 # Number of rows
nCol = 50 # Number of columns
p = range(1, 10, 1) # range of values of p (proportion of NLC)
h = range(0, 11, 1) # range of values of h (spatial autocorrelation)
r = range(10) # number of replicates per landscape type
w = [3,5,9,15] # window sizes for moving window analysis
ex_keys = [{'a': 0.0000001}, {'a':1}, {'a':5}, {'a':50}, {'a':-1}, {'a':-5}, {'a':-50}]

# create the parameters to create surfaces for based on how they were set up above
param_set = pd.DataFrame(expandgrid(p, h, r, w, ex_keys))
param_set = param_set.rename(index=str, columns = {'Var1':'p', 'Var2':'h', 'Var3':'r', 'Var4':'w', 'Var5':'ex_keys'})
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