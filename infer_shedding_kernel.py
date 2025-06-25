import pandas as pd
import pickle
import json
import requests
from scipy import signal
from datetime import date,timedelta
import yaml
import copy
import numpy as np

import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

start = pd.to_datetime('2023-01-01')
eval_dates = [start + datetime.timedelta(days=7*j) for j in range(0,52)]

cases = pd.read_csv('data/covid19cases_test.csv',index_col=0)
cases = cases[cases['area'] =='San Diego']
cases = cases[(cases.index>='2020-06-01') & (cases.index<='2022-07-01')]
cases.index = pd.to_datetime(cases.index)

cases['cases'] = cases['cases'].rolling(window=7, center=True, min_periods=0).mean()
cases['positivity'] = cases['cases'].rolling(window=7, center=True, min_periods=0).mean()/cases['total_tests'].rolling(window=7, center=True, min_periods=0).mean()

ww_dict = {'Point Loma':['data/PointLoma_sewage_seqs.csv','data/PointLoma_sewage_qPCR.csv']}
for site, files in zip(ww_dict.keys(),ww_dict.values()):
    df = pd.read_csv(f'{files[0]}')

    df['Date'] = pd.to_datetime(df['Date'])
    df.columns = [dfc.split(' (')[0] for dfc in df.columns]
    df =df.set_index('Date')
    df = df[df.index>='2020-06-01']
    df = df[df.index<='2022-07-01']
    df = df.dropna(axis = 0, how = 'all')
    df = df.fillna(0)
    df = df/100.

    df = df.drop(columns=['Other'])
    df = df[df.columns[df.sum(axis=0) > 0.01]]

    cdf = pd.read_csv(f'{files[1]}')
    cdf['Sample_Date'] = pd.to_datetime(cdf['Sample_Date'])
    cdf =cdf.set_index('Sample_Date')
    sharedInds = np.sort(list(set(cdf.index) & set(df.index)))
    cdf = cdf.loc[sharedInds]
    df = df.loc[sharedInds]
    df = df[~df.index.duplicated(keep='last')]
    scaleddf = df.mul(cdf['Mean viral gene copies/L'],axis=0)

    cdf = cdf.resample('D').asfreq()
    cdf = cdf.rolling(window=7, center=True, min_periods=0).mean()


cdf = pd.concat([cdf,cases],axis=1)
cdf = cdf.dropna(how='any')

N = 30
F = 5
# first, let's learn the shedding kernel. 
X = np.array([cdf['Mean viral gene copies/L'].values[(j-F):(j+N-F)] for j in range(F,(cdf.shape[0]-N+F))])/cdf['Mean viral gene copies/L'].mean()
#add coeficient term
Y= cdf['cases'].values[F:(len(cdf['cases'])-N+F)]

import cvxpy as cp

fig,ax = plt.subplots()
ax.plot(cdf.index,cdf['Mean viral gene copies/L'],color='cornflowerblue')
ax.set_ylabel('Mean viral gene copies/L')
ax2 = ax.twinx()
ax2.plot(cdf.index,cdf['cases'],color='black')
ax2.set_ylabel('Clinical cases')
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
fig.tight_layout()
plt.savefig('plots/epi_curves_plot.pdf',bbox_inches='tight')

# least squares problem
x = cp.Variable(N)
cost = cp.norm(X @ x - Y,2)
constraints = [x >= 0]
prob = cp.Problem(cp.Minimize(cost),constraints)
prob.solve(verbose=True)#,solver=cp.CLARABEL)

Y2= cdf['positivity'].values[0:(len(cdf['positivity'])-N)]

#switch to logistic function
beta = cp.Variable(N,nonneg=True)#+1)
log_likelihood = cp.sum(
    cp.multiply(Y2, X @ beta) - cp.logistic(X @ beta)
)
# constraints = [beta >= 0.001, beta<=10]
prob2 = cp.Problem(cp.Maximize(log_likelihood))#,constraints)
prob2.solve(verbose=True)#,solver=cp.CLARABEL)

fig,ax = plt.subplots()
ax.plot(range(-F,len(x.value)-F),x.value) # just plotting the curve learned via the case counts for now. 
ax.set_xlabel('Date relative to detection')
ax.set_ylabel('Mean viral load copies/L')
plt.savefig('plots/learned_shedding_kernel.pdf')
