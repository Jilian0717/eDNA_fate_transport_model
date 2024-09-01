import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
import pickle, os, cmocean
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.dates as mdates
import seaborn as sns
import scipy

#% 3 models
import pandas as pd
fn = pd.read_csv('eDNA_distance_particle_rpk.csv')
eDNA = fn['ModifiedSiteMean'] # original eDNA + 0.1 copiesL-1
eDNA = np.log(eDNA)
dist = fn['Distance']
nowind = fn['NoWind'] # no wind particle percentages
realistic = fn['Realistic']
detection = fn['Detected']

fig = plt.figure(figsize=(16,8))
plt.style.use('ggplot')
fz = 17; sz=100; alpha=0.2
ax = fig.add_subplot(231)
plt.scatter(dist, detection, s=sz, c='k', alpha=alpha)
plt.scatter(dist, detection, s=sz,edgecolors='k', facecolors="none", linewidth=1)
plt.xlabel('Distance from source (km)', fontsize=fz)
plt.ylabel('Probability of Detection', fontsize=fz)
ax.tick_params(labelsize=fz)
plt.title('(a) Distance-only model', x=0.3, fontsize=18)

ax = fig.add_subplot(232)
plt.scatter(nowind, detection, s=sz, c='k', alpha=alpha)
plt.scatter(nowind, detection, s=sz, edgecolors='k', facecolors="none", linewidth=1)
plt.xlabel('Simulated particles (no wind, %)', fontsize=fz)
ax.tick_params(labelsize=fz)
plt.title('(b) No-wind model', x=0.22, fontsize=18)

ax = fig.add_subplot(233)
plt.scatter(realistic, detection, s=sz, c='k', alpha=alpha)
plt.scatter(realistic, detection, s=sz, edgecolors='k', facecolors="none", linewidth=1)
plt.xlabel('Simulated particles (%)', fontsize=fz)
ax.tick_params(labelsize=fz)
# load logistic fitting line
fn1 = pd.read_csv('logistic_model_fittedvalues.csv')
y_low_95ci = fn1.low_95ci.values
y_high_95ci = fn1.high_95ci.values
xfit = fn1['par %'].values; xfit_sort=np.sort(xfit)
yfit = fn1['mean_est'].values
sort_idx = np.argsort(xfit)
plt.plot(xfit_sort, yfit[sort_idx], c='k', lw=1.7)
plt.fill_between(xfit_sort, y_low_95ci[sort_idx], y_high_95ci[sort_idx], 
                 color='k', alpha=0.15, edgecolor='none')
plt.title('(c) eDNA fate and transport model', x=0.45, fontsize=18)

ax = fig.add_subplot(234)
plt.scatter(dist, eDNA, s=sz, c='k', alpha=alpha)
plt.scatter(dist, eDNA, s=sz, edgecolors='k', facecolors="none", linewidth=1)
plt.xlabel('Distance from source (km)', fontsize=fz)
plt.ylabel('log(copies + 0.1 L$^{-1}$)', fontsize=fz)
ax.tick_params(labelsize=fz)
plt.title('(d) Distance-only model', x=0.3, fontsize=18)
add_linear_regression = False
if add_linear_regression:
    x = dist
    y = eDNA
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    ax.text(0.75,0.19, 'R$^2$='+"%.2f" % (r_value**2), transform=ax.transAxes, fontsize=15)
    ax.text(0.75,0.10, 'p='+"%.4f" % (p_value), transform=ax.transAxes, fontsize=15)

ax = fig.add_subplot(235)
plt.scatter(nowind, eDNA, s=sz, c='k', alpha=alpha)
plt.scatter(nowind, eDNA, s=sz, edgecolors='k', facecolors="none", linewidth=1)
plt.xlabel('Simulate particles (no wind, %)', fontsize=fz)
ax.tick_params(labelsize=fz)
plt.title('(e) No-wind model', x=0.22, fontsize=18)
add_linear_regression = False
if add_linear_regression:
    x = nowind
    y = eDNA
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    ax.text(0.75,0.19, 'R$^2$='+"%.2f" % (r_value**2), transform=ax.transAxes, fontsize=15)
    ax.text(0.75,0.10, 'p='+"%.4f" % (p_value), transform=ax.transAxes, fontsize=15)


ax = fig.add_subplot(236)
plt.scatter(realistic, eDNA, s=sz, c='k', alpha=alpha)
plt.scatter(realistic, eDNA, s=sz, edgecolors='k', facecolors="none", linewidth=1)
plt.xlabel('Simulated particles (%)', fontsize=fz)
ax.tick_params(labelsize=fz)
# load 95% band
fn = pd.read_csv('linear_model_fittedvalues.csv')
yfit = fn.mean_est.values
y_high_95ci = fn.high_95ci.values
y_low_95ci = fn.low_95ci.values
xfit = fn['par %'].values
xfit_sort=np.sort(xfit)
sort_idx = np.argsort(xfit)
plt.plot(xfit,yfit,'k',lw=2)
plt.fill_between(xfit_sort, y_low_95ci[sort_idx], y_high_95ci[sort_idx], 
                 color='k', alpha=0.15, edgecolor='none')
plt.title('(f) eDNA fate and transport model', x=0.45, fontsize=18)
x = realistic
y = eDNA
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
plt.tight_layout()
#plt.savefig('three-models', dpi=300, bbox_inches='tight')
