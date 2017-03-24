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

if not os.path.exists('two_step_binary/results'):
	os.makedirs('two_step_binary/results')

p_vals = h_vals = np.arange(0.1, 0.9, 0.1)
w1_vals = w2_vals = [3, 9, 15]
f1_vals = f2_vals = ["linear", "exp", "negexp"]

res = pd.DataFrame()

# horrible horrible for loop. 
for p in p_vals:
    for h in h_vals:
        for w1 in w1_vals:
            for w2 in w2_vals:
                for f1 in f1_vals:
                    for f2 in f2_vals:
                        out = nlm.two_step_binary(50, p, h, w1, w2, f1, f2)
                        res = res.append(out, ignore_index=True)                           
                        
unique_filename = uuid.uuid4()
res.to_csv('two_step_binary/results/output'+str(unique_filename)+'.csv', index=False)




