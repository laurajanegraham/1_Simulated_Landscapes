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

if not os.path.exists('farmland_birds/results'):
	os.makedirs('farmland_birds/results')

h_vals = np.arange(0, 1, 0.1)
w_vals = [3, 9, 15]
npp_vals = np.arange(0.1, 1, 0.1)

res = pd.DataFrame()

# horrible horrible for loop. 
for h in h_vals:
    for w in w_vals:
        for npp in npp_vals:
            out = nlm.farmland_birds(50, h, w, w, npp)
            res = res.append(out, ignore_index=True)                           
                        
unique_filename = uuid.uuid4()
res.to_csv('farmland_birds/results/output'+str(unique_filename)+'.csv', index=False)

