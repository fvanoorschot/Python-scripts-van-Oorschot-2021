# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:22:00 2021

@author: fransjevanoorschot

TEST SR_CALCULATION FUNCTIONS

Sr calculation with memory method as described in paper

function 1: interception calculation and initial storage deficit calculation
function 2: iterative approach (10 iterations) for calculating Et and Sd
function 3: Sr calculation based on return periods

"""

#%% LOAD PACKAGES
import numpy as np
import pandas as pd
from scipy.stats import gumbel_r
import os,glob
import sr_calculation as sr

#%% test
date_start = '03-01'
date_end = '02-28'
year_start = 1973
year_end = 1983
Si_0 = 0  
Si_max = 2.5 
it = 3
#csv with catchment daily P Q and Ep
filename = 'C:/Users/fransjevanoors/surfdrive/Fransje/Thesis/Output_data/All_data/2706_WA/daily/Abercrombie_7310_D.csv'

# run function 1
A = sr.sd_initial(filename, Si_0, Si_max,year_start,year_end,date_start,date_end)

# run function 2
it_out = sr.sd_iterations(A,date_start,date_end,year_start,year_end,it)

# run function 3
T=[2,5,10,15,20,30,40]
Sd = it_out[3]
rp_out = sr.sr_return_periods(T,Sd,date_start,date_end,year_start,year_end,it)

