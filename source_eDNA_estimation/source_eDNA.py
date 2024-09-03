import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.dates as mdates

# plot source concentration # passive particle 100k, r=100m, z=2m
fn = 'estimate_source_eDNA_new.xlsx'
tmp = pd.read_excel(fn, header=0)
t_r100m = tmp['shedding time'].values
C0_r100m = tmp['C0-r100m'].values
C0_r40m = tmp['C0-r40m'].values
C0_r60m = tmp['C0-r60m'].values
C0_r80m = tmp['C0-r80m'].values
C0_r120m = tmp['C0-r120m'].values
C0_r140m = tmp['C0-r140m'].values
C0_r160m = tmp['C0-r160m'].values
C0_r180m = tmp['C0-r180m'].values
C0_r200m = tmp['C0-r200m'].values
sample_time = tmp['Sampling time (UTC)'].values
eDNA_8site = tmp['eDNA_measured (copies L-1)'].values  # 8 solid stations with positive DNA and large particles
eDNA_8site_std = tmp['std'].values
stations8 = tmp['Site'].values
travel_t = tmp['particle travel time (hour) - r100m']
travel_d = tmp['travel_d(km)']
#% source eDNA 
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
mpl.rcParams['axes.linewidth'] = 1
for i in range(len(stations8)):
    test_r = True
    if test_r:
        plt.plot(t_r100m[i], C0_r40m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.1)
        plt.plot(t_r100m[i], C0_r60m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.2)
        plt.plot(t_r100m[i], C0_r80m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.3)
        plt.plot(t_r100m[i], C0_r100m[i], 'o', markeredgecolor='k', markerfacecolor='b',ms=8, markeredgewidth=0.5, alpha=0.9)
        plt.plot(t_r100m[i], C0_r120m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.5)
        plt.plot(t_r100m[i], C0_r140m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.6)
        plt.plot(t_r100m[i], C0_r160m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.7)
        plt.plot(t_r100m[i], C0_r180m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.8)
        plt.plot(t_r100m[i], C0_r200m[i], 'o', markeredgecolor='k', markerfacecolor='g',ms=8, markeredgewidth=0.5, alpha=0.9)
    else:
        plt.plot(t_r100m[i], C0_r100m[i], '*', markeredgecolor='k', markerfacecolor='g', 
             ms=15, markeredgewidth=0.5, alpha=0.7)
    if stations8[i][-2:]=='15' or stations8[i][-2:]=='39' :
        plt.text(t_r100m[i]-np.timedelta64(1,'h'), C0_r100m[i], 'B'+stations8[i][-2:], fontsize=12, c='k')
    elif stations8[i][-2:]=='41' :
        plt.text(t_r100m[i]+np.timedelta64(10,'m'), C0_r100m[i]*0.4, 'B'+stations8[i][-2:], fontsize=12, c='k')
    else:
        plt.text(t_r100m[i], C0_r100m[i]*1.01, 'B'+stations8[i][-2:], fontsize=12, c='k')
    # measured
    errorbar=True
    if errorbar:
        err = plt.errorbar(sample_time[i], eDNA_8site[i], fmt="o", markersize=8,
                  yerr = eDNA_8site_std[i], color='brown', lw=1.2,alpha=0.9,capsize=4)
        err[-1][0].set_linestyle('--')
    else:
        plt.plot(sample_time[i], eDNA_8site[i], '*', markeredgecolor='brown', markerfacecolor='brown', 
             ms=15, markeredgewidth=0.5, label='measured', alpha=0.9)
    if stations8[i][-2:]=='15':
        plt.text(sample_time[i]-np.timedelta64(30,'m'), eDNA_8site[i]*0.1, 'B'+stations8[i][-2:], fontsize=12, c='k')
    elif stations8[i][-2:]=='16':
        plt.text(sample_time[i], eDNA_8site[i]*1.8, 'B'+stations8[i][-2:], fontsize=12, c='k')
    elif stations8[i][-2:]=='14':
        plt.text(sample_time[i]+np.timedelta64(10,'m'), eDNA_8site[i]*1.2, 'B'+stations8[i][-2:], fontsize=12, c='k')
  
    elif stations8[i][-2:]=='40':
        plt.text(sample_time[i]-np.timedelta64(1,'h'), eDNA_8site[i]*1, 'B'+stations8[i][-2:], fontsize=12, c='k')
    elif stations8[i][-2:]=='39':
        plt.text(sample_time[i]-np.timedelta64(1,'h'), eDNA_8site[i]*0.4, 'B'+stations8[i][-2:], fontsize=12, c='k')
    elif stations8[i][-2:]=='41':
        plt.text(sample_time[i]+np.timedelta64(10,'m'), eDNA_8site[i]*0.8, 'B'+stations8[i][-2:], fontsize=12, c='k')
  
    else:
        plt.text(sample_time[i], eDNA_8site[i], 'B'+stations8[i][-2:], fontsize=12, c='k')
    plt.yscale('log')
    
# add a box for sampling campaign
rect = patches.Rectangle((datetime(2023,6,5,18,40,0), 0), timedelta(seconds=3600*4.6), 
                         1e3, linewidth=1, edgecolor='brown', facecolor='brown', alpha=0.05)
ax.add_patch(rect)
# add a box for shedding
if test_r:
    rect = patches.Rectangle((datetime(2023,6,5,0,0,0), 6*1e3), timedelta(seconds=3600*12), 
                         2*1e7, linewidth=1, edgecolor='g', facecolor='g', alpha=0.05)
else:
    rect = patches.Rectangle((datetime(2023,6,5,0,0,0), 2e4), timedelta(seconds=3600*12), 
                         1e7, linewidth=1, edgecolor='g', facecolor='g', alpha=0.05)
ax.add_patch(rect)
plt.text(datetime(2023,6,5,4,0,0), 6e6, 'Source eDNA', fontsize=14, color='g', alpha=0.8, fontweight='bold')
plt.text(datetime(2023,6,5,18,0,0), 3*1e3, 'Measured eDNA', fontsize=14, c='brown', fontweight='bold')
plt.xticks(rotation=0, fontsize=16 ) #,  labels=['06:00','12:00','18:00','00:00','06:00','12:00','18:00','00:00'])
plt.yticks(fontsize=18)
plt.xlabel('June 5-6, 2023', fontsize=20)
plt.ylabel('eDNA copies L$^{-1}$', fontsize=20)
myFmt = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(myFmt)
plt.tight_layout()
#plt.savefig('initial_eDNA_from_particle_v1',dpi=300, bbox_inches='tight')
