# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:16:00 2017

@author: lg1u16
"""
import numpy as np
import pandas as pd
import os
import nlmfunctions as nlm
import uuid

if not os.path.exists('two_step_binary_cont/results'):
	os.makedirs('two_step_binary_cont/results')

p_vals = h_vals = np.arange(0.1, 0.9, 0.1)
h2_vals = [0, 1]
w1_vals = w2_vals = w3_vals = [3, 9, 15]

res = pd.DataFrame()

# horrible horrible for loop. 
for p in p_vals:
    for h1 in h_vals:
        for w1 in w1_vals:
            for w2 in w2_vals:
                for w3 in w3_vals:
                    for h2 in h2_vals:
                        out = nlm.two_step_binary_cont(50, p, h1, h2, w1, w2, w3)
                        res = res.append(out, ignore_index=True)                           
                        
unique_filename = uuid.uuid4()
res.to_csv('two_step_binary_cont/results/output'+str(unique_filename)+'.csv', index=False)

