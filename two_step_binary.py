# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:16:00 2017

@author: lg1u16
"""
import pandas as pd
import os

import sys


import nlmfunctions as nlm


# gets the number of the job array - use this to select the line from the params file
job_no = sys.argv[1]
job_no = int(job_no)

params = pd.read_csv("two_step_binary/params.csv")

params = params.ix[job_no]

out = nlm.two_step_binary(50, params['p'], params['h'], params['w1'], params['w2'], params['f1'], params['f2'])

if not os.path.exists('two_step_binary/results'):
	os.makedirs('two_step_binary/results')

out.to_csv('two_step_binary/results/output'+str(job_no)+'.csv', index=False)

