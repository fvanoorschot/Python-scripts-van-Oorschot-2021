# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 14:53:09 2020

@author: fransjevanoors

SIGNIFICANCE CALCULATIONS FUNCTION

Approach:
1. test significance of correlation improvement
2. H0: each couple of datapoints (at every timestep) come from the same population; there is no difference between model1 and model2
3. Randomly select for each timestep either model1 or model2 and generate 2x1000 random samples
4. Calculate correlation coefficients for random samples with observations
5. Calculate correlation differences between randomly selected random samples
6. Distribution of correlation differences -> 5% and 10% tests (one tailed because only improvement is interesting)

"""
#%% LOAD PACKAGES
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import norm

#%% FUNCTION
def statistics_correlation_difference(obs, data_1, data_2, alpha, reps):
    
    # correlation for data1 and data2 with observations
    r_1 = pearsonr(obs,data_1)[0]
    r_2 = pearsonr(obs,data_2)[0]
    r_diff = r_1 - r_2
    
    n = len(obs)
    
    # make 1000 artificial datasets
    x = np.zeros((reps,n))
    y = np.zeros((reps,n))
    for j in range(reps):
        for i in range(n):
            array = [data_1.iloc[i],data_2.iloc[i]]
            x[j,i] = np.random.choice(array,size=1)
            y[j,i] = np.random.choice(array,size=1)
    
    # calculate correlation
    r_x = np.zeros(reps)
    r_y = np.zeros(reps)
    for j in range(reps):
        r_x[j] = pearsonr(x[j,:],obs)[0]
        r_y[j] = pearsonr(y[j,:],obs)[0]

    # randomly generate correlation difference sample
    r_diff_sample = np.random.choice(r_x,size=reps) - np.random.choice(r_y,size=reps)
    
    # plot distribution of r
    plt.figure(figsize=(10,6))
    ci = np.percentile(r_diff_sample,alpha)    
    plt.hist(r_diff_sample,reps,histtype='step',label='random sample rdiff')
    plt.axvline(ci,color='k',linestyle='-.',label=str(alpha)+'%')
    plt.axvline(r_diff,color='r',label='r1 - r2')
    plt.xlabel('Correlation difference')
    plt.legend()
    
    return(r_diff,ci)
    

#%% RUN FUNCTION - example
folder = 'C:/Users/fransjevanoors/surfdrive/Fransje/HTESSEL/HTESSEL_output_aug2020/MAL3L4/MM_dataframes/'
# a is dataframe with monthly observed and modeled (CTR and MD models) discharge Q for a specific catchment
a = pd.read_csv(str(folder)+'Q_M_IA_EA.csv',index_col=0) 

obs = a.obs #Observations
data_1 = a.base #CTR model
data_2 = a.MD #MD model
reps = 1000 #repetitions
alpha = 5 #5 percent confidence interval

statistics_correlation_difference(obs, data_1, data_2, alpha, reps)

