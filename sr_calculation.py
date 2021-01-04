# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 09:55:19 2020

@author: fransjevanoorschot

Sr calculation with memory method as described in paper

function 1: interception calculation and initial storage deficit calculation
function 2: iterative approach (10 iterations) for calculating Et and Sd
function 3: Sr calculation based on return periods

"""

#%% LOAD PACKAGES
import numpy as np
import pandas as pd
from scipy.stats import gumbel_r

#%% function 1: interception calculation and intial storage deficits

# INPUT
# filename: csv file with daily catchment values for P, Ep, Q
# Si_0: initial interception storage = 0
# Si_max: maximum interception storage = 2.5mm
# date_start, date_end: start and end 'month-day' of time-series (depending on hydro-year)
# year_start, year_end: start and end year of time-series

# OUTPUT
# catchment: pandas dataframe with daily catchment values for P, Ep, Q, Pe, Et and Sd (based on initial Et estimate)

def sd_initial(filename, Si_0, Si_max,year_start,year_end,date_start,date_end):
    
    #read csv file for catchment of interest
    catchment = pd.read_csv(filename, sep=',', skiprows=0, index_col=0, skipinitialspace=True)
    catchment = catchment.loc[str(year_start)+'-'+str(date_start):str(year_end)+'-'+str(date_end)]
    catchment.index = pd.to_datetime(catchment.index)
    # add columns for interception storage calculation
    catchment['Si_1'] = np.nan
    catchment['Pe'] = np.nan
    catchment['Si_2'] = np.nan
    catchment['Ei'] = np.nan
    catchment['Si_3'] = np.nan
    catchment['Et'] = np.nan
    catchment['Sd'] = np.nan
    
    #calculate interception storage and effective precipitation for all timesteps
    for i in range(1,np.size(catchment.index)):
        catchment.Si_1[0] = catchment.P[0] + Si_0
        catchment.Pe[0] = max(0,catchment.Si_1[0]-Si_max)
        catchment.Si_2[0] = catchment.Si_1[0] - catchment.Pe[0]
        catchment.Ei[0] = min(catchment.Si_2[0],catchment.Ep[0])
        catchment.Si_3[0] = catchment.Si_2[0] - catchment.Ei[0]
    
        catchment.Si_1[i] = catchment.P[i] + catchment.Si_3[i-1]
        catchment.Pe[i] = max(0,catchment.Si_1[i]-Si_max)
        catchment.Si_2[i] = catchment.Si_1[i] - catchment.Pe[i]
        catchment.Ei[i] = min(catchment.Si_2[i],catchment.Ep[i])
        catchment.Si_3[i] = catchment.Si_2[i] - catchment.Ei[i]

    #water balance Et calculation (Et = Pe-Q)
    Pe_sum = np.sum(catchment.Pe)
    Ep_sum = np.sum(catchment.Ep)
    Q_sum = np.sum(catchment.Q)
    Et_sum = Pe_sum - Q_sum
    
    #calculate daily Et (Ep(daily)*(Et_sum/Ep_sum)) and Sd
    for i in range(1,np.size(catchment.index)):
        catchment.Et[0] = catchment.Ep[0]/Ep_sum * Et_sum
        catchment.Sd[0] = catchment.Pe[0] - catchment.Et[0]
    
        catchment.Et[i] = catchment.Ep[i]/Ep_sum * Et_sum
        catchment.Sd[i] = min(0,catchment.Sd[i-1]+catchment.Pe[i]-catchment.Et[i])
    
    return catchment

#%% function 2: iterative approach (10 iterations) for calculating Et and Sd

# INPUT
# A: pandas dataframe with daily data for P, Ep, Q and initial Et and Sd (output from sd_initial function)
# date_start, date_end: start and end 'month-day' of time-series (depending on hydro-year)
# year_start, year_end: start and end year of time-series
# it: amount of iterations
    
# OUTPUT
# output_a: annual values for Q, Pe, Ep and Et
# c: inter-annual variation coefficient for all years and all iterations
# Et: calculated daily transpiration for all iterations
# Sd: calculated daily storage deficits for all iterations
# Et_a_it: annual transpiration for all yeras and all iterations
# WB: water balance for all iterations
    
def sd_iterations(A,date_start,date_end,year_start,year_end,it):
    
    # make dataframes where iteration results are stored
    Sd = pd.DataFrame(index=pd.date_range(start=str(year_start)+'-'+str(date_start), end=str(year_end)+'-'+str(date_end)),columns=['Sd0'])
    Et = pd.DataFrame(index=pd.date_range(start=str(year_start)+'-'+str(date_start), end=str(year_end)+'-'+str(date_end)),columns=['Et0'])
    Et.Et0 = A.Et
    Sd.Sd0 = A.Sd
    
    # arrays with daily Ep and Pe values
    A_Ep = np.zeros(len(A.Ep))
    A_Ep = (A.Ep.loc[str(year_start)+'-'+str(date_start):str(year_end)+'-'+str(date_end)])
    A_Pe = np.zeros(len(A.Ep))
    A_Pe = (A.Pe.loc[str(year_start)+'-'+str(date_start):str(year_end)+'-'+str(date_end)])
    
    #empty arrays for annual sums of Q, Pe, Ep and Et
    total_years = year_end-year_start
    Q_a = np.zeros(total_years)
    Pe_a = np.zeros(total_years)
    Ep_a = np.zeros(total_years)
    Et_a = np.zeros(total_years)
    
    years = range(year_start,year_end+1,1)
    for i in range(0,total_years,1):
        Q_a[i]= (A.Q.loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)].sum())
        Pe_a[i]= (A.Pe.loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)].sum())
        Ep_a[i]= (A.Ep.loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)].sum())
        Et_a[i]= (A.Et.loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)].sum())
    
    #store annual sums in dataframe
    output_a = pd.DataFrame(Q_a, columns=['Q_a'])
    output_a['Pe_a'] = Pe_a
    output_a['Ep_a'] = Ep_a
    output_a['Et_a'] = Et_a
    output_a.index = pd.date_range(start=str(year_start)+'-'+str(date_start),periods=total_years,freq='12M')
    
    #dataframe with annual values for Et for all iterations
    Et_a_it = pd.DataFrame(columns=['Et0'])
    Et_a_it.Et0 = Et_a #first iteration is initial Et
    
    #dataframe with inter-annual variability coefficient c for all iterations
    c = pd.DataFrame(index=np.arange(0,total_years,1),columns=['c0'])
    c.c0 = (sum(Et_a))/(sum(Ep_a)) #initial c
    c0 = (sum(Et_a))/(sum(Ep_a))

    #ITERATION 1
    k=1
    dSdt = np.zeros(total_years)
    #add columns to dataframes
    Sd['Sd'+str(k)] = np.nan
    c['c'+str(k)] = np.nan
    Et['Et'+str(k)] = np.nan
    Et_a_it['Et'+str(k)] = np.nan
    #calculate change in storage dSdt for each year
    for i in range(total_years):
        dSdt[i] = (Sd['Sd'+str(k-1)].loc[str(years[i])+'-'+str(date_start)]-Sd['Sd'+str(k-1)].loc[str(years[i+1])+'-'+str(date_end)])
    
    #annual water balance and annual value for c
    Et_a1 = Pe_a - Q_a - dSdt
    c['c'+str(k)] = Et_a1/Ep_a
    
    #set constraints on c, should be within +-25% of c0
    for i in range(np.size(c['c'+str(k)])):
        if c['c'+str(k)][i]>= c0*1.25:
            c['c'+str(k)][i]=c0*1.25
        if c['c'+str(k)][i] <= c0*0.75:
            c['c'+str(k)][i]=c0*0.75
    
    #calculate Et with c
    for j in range(total_years):
        Et['Et'+str(k)].loc[str(years[j])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)]=c['c'+str(k)].loc[j]*A.Ep.loc[str(years[j])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)]
        Et_a_it['Et'+str(k)].loc[j] =c['c'+str(k)].loc[j]*Ep_a[j] 

    #set constraints to WB, in case more than 2% deviation -> distribute excess/lacking water over the time-series
    WB = sum(Pe_a)-sum(Q_a)-sum(Et['Et'+str(k)])
    if np.abs(WB) > 0.02*sum(Pe_a):
        Et['Et'+str(k)] = Et['Et'+str(k)] * ((WB+sum(Et['Et'+str(k)]))/sum(Et['Et'+str(k)])) 
        for i in range(0,total_years,1):
            Et_a_it['Et'+str(k)].loc[i]=Et['Et'+str(k)].loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)].sum()
    
    #calculate storage deficits for full timeseries        
    for j in range(1,np.size(Et.index)):
        Sd['Sd'+str(k)][0] = A_Pe[0]-Et['Et'+str(k)][0]
        Sd['Sd'+str(k)][j] = min(0,Sd['Sd'+str(k)][j-1]+A_Pe[j]-Et['Et'+str(k)][j])
 
    #ITERATION 2 TO end(it), TAKE AVERAGE OF PREVIOUS 2 ITERATIONS FOR Sd
    for k in range(2,it,1):
        #add columns to dataframes
        Sd['Sd'+str(k)] = np.nan
        c['c'+str(k)] = np.nan
        Et['Et'+str(k)] = np.nan
        Et_a_it['Et'+str(k)] = np.nan
        
        #calculate change in storage dSdt for each year
        dSdt = np.zeros(total_years)
        for i in range(total_years):
            Sd0 = (Sd['Sd'+str(k-2)].loc[str(years[i])+'-'+str(date_start)]+Sd['Sd'+str(k-1)].loc[str(years[i])+'-'+str(date_start)])/2
            Sd1 = (Sd['Sd'+str(k-2)].loc[str(years[i+1])+'-'+str(date_end)]+Sd['Sd'+str(k-1)].loc[str(years[i+1])+'-'+str(date_end)])/2
            dSdt[i] = Sd0-Sd1 
        
        #annual water balance and annual value for c
        Et_a1 = Pe_a - Q_a - dSdt
        c['c'+str(k)] = Et_a1/Ep_a
        
        #set constraints on c, should be within +-25% of c0
        for i in range(np.size(c['c'+str(k)])):
            if c['c'+str(k)][i]>= c0*1.25:
                c['c'+str(k)][i]=c0*1.25
            if c['c'+str(k)][i] <= c0*0.75:
                c['c'+str(k)][i]=c0*0.75
        
        #calculate Et with c
        for j in range(total_years):
            Et['Et'+str(k)].loc[str(years[j])+'-'+str(date_start):str(years[j+1])+'-'+str(date_end)]=c['c'+str(k)].loc[j]*A.Ep.loc[str(years[j])+'-'+str(date_start):str(years[j+1])+'-'+str(date_end)]
            Et_a_it['Et'+str(k)].loc[j] =c['c'+str(k)].loc[j]*Ep_a[j] 
        
        #set constraints to WB, in case more than 2% deviation -> distribute excess/lacking water over the time-series     
        WB = sum(Pe_a)-sum(Q_a)-sum(Et['Et'+str(k)])
        if np.abs(WB) > 0.02*sum(Pe_a):
            Et['Et'+str(k)] = Et['Et'+str(k)] * ((WB+sum(Et['Et'+str(k)]))/sum(Et['Et'+str(k)])) 
            for i in range(0,total_years,1):
                Et_a_it['Et'+str(k)].loc[i]=Et['Et'+str(k)].loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end)].sum()
        
        #calculate storage deficits for full timeseries        
        for j in range(1,np.size(Et.index)):
            Sd['Sd'+str(k)][0] = A_Pe[0]-Et['Et'+str(k)][0]
            Sd['Sd'+str(k)][j] = min(0,Sd['Sd'+str(k)][j-1]+A_Pe[j]-Et['Et'+str(k)][j])
            
    #store values for annual water balances        
    WB = []
    for i in range(1,it,1):
        WB.append(sum(Pe_a)-sum(Q_a)-sum(Et_a_it['Et'+str(i-1)]))
                
    return(output_a,c,Et,Sd,Et_a_it,WB)

#%% function 3: Sr calculation based on return periods

# INPUT
# T: array of return periods of interest T=[2,5,10,15,20,30,40]
# Sd: dataframe of Sd calculated in sd_iterations function
# date_start, date_end: start and end 'month-day' of time-series (depending on hydro-year)
# year_start, year_end: start and end year of time-series
# it: amount of iterations
    
# OUTPUT
# Sd_T: storage deficits corresponding with return periods T
    

def sr_return_periods(T,Sd,date_start,date_end,year_start,year_end,it):

    for j in range(len(T)):
        Sd = Sd*-1
        total_years = year_end-year_start
        years = range(year_start,year_end+1,1)
        
        # calculate annual max Sd
        Sd_max=[]
        for i in range(0,total_years,1):
            Sd_max.append(max(Sd.loc[str(years[i])+'-'+str(date_start):str(years[i+1])+'-'+str(date_end),'Sd'+str(it-1)]))
                   
        # gumbel function
        def gumbel_r_mom(x):
            scale = np.sqrt(6)/np.pi * np.std(x)
            loc = np.mean(x) - np.euler_gamma*scale
            return loc, scale    
        
        loc1, scale1 = gumbel_r_mom(Sd_max)
        
        # find Sd value corresponding with return period
        Sd_T = []
        for i in np.arange(0,len(T),1):
            p = 1-(1/T[i])
            y = -np.log(-np.log(p))
            x = scale1 * y + loc1
            Sd_T.append(x)
         
        return(Sd_T)   

