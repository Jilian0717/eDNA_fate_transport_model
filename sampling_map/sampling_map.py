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

# grid file
grd='grid.nc'
dsg = xr.open_dataset(grd)
latr = dsg['lat_rho'].values
lonr = dsg['lon_rho'].values
h = dsg['h'].values
maskr = dsg['mask_rho'].values
hm = np.ma.masked_where(maskr==1, h)
hm2 = np.ma.masked_where(maskr==0, h)

# sample info by sampling order-16, 17, 15, 14, 35, 36, 37, 38, 39, 40, 41, 42, 43, 09, 11
stations = dict()
stations['BANG16'] = {'lat': 47.7743587, 'lon': -122.7072554, 'sample_time0': '1051-1057', 
                      'conc': np.nan, 'sd':np.nan,  # copies/L, a place-holder
                      'model_time': datetime(2023,6,5,19,0,0), 'sample_t_local': datetime(2023,6,5,10,51,0)}
stations['BANG17'] = {'lat': 47.77476498, 'lon': -122.7139617, 'sample_time0': '1103-1109', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,19,0,0), 'sample_t_local': datetime(2023,6,5,11,3,0)}
stations['BANG15'] = {'lat': 47.77071505, 'lon': -122.7132016, 'sample_time0': '1118-1129', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,19,30,0), 'sample_t_local': datetime(2023,6,5,11,18,0)}
stations['BANG14'] = {'lat': 47.77053826, 'lon': -122.7181893, 'sample_time0': '1138-1144', 
                      'conc': np.nan, 'sd':np.nan, 
                      'model_time': datetime(2023,6,5,19,30,0), 'sample_t_local': datetime(2023,6,5,11,38,0)}
stations['BANG35'] = {'lat': 47.76205, 'lon': -122.727267, 'sample_time0': '1155-1202', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,20,0,0), 'sample_t_local': datetime(2023,6,5,11,55,0)}
stations['BANG36'] = {'lat': 47.758383, 'lon': -122.732183, 'sample_time0': '1209-1219', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,20,0,0), 'sample_t_local': datetime(2023,6,5,12,9,0)}
stations['BANG37'] = {'lat': 47.755, 'lon': -122.736083, 'sample_time0': '1225-1235', 
                      'conc': np.nan,'sd':np.nan,
                      'model_time': datetime(2023,6,5,20,30,0), 'sample_t_local': datetime(2023,6,5,12,25,0)}
stations['BANG38'] = {'lat': 47.747767, 'lon': -122.741083, 'sample_time0': '1244-1257', 
                      'conc': np.nan,'sd':np.nan,
                      'model_time': datetime(2023,6,5,20,30,0), 'sample_t_local': datetime(2023,6,5,12,44,0)}
stations['BANG39'] = {'lat': 47.742333, 'lon': -122.742417, 'sample_time0': '1305-1318', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,21,0,0), 'sample_t_local': datetime(2023,6,5,13,5,0)}
stations['BANG40'] = {'lat': 47.738267, 'lon': -122.74495, 'sample_time0': '1323-1338', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,21,30,0), 'sample_t_local': datetime(2023,6,5,13,23,0)}
stations['BANG41'] = {'lat': 47.7297, 'lon': -122.749283, 'sample_time0': '1355-1405', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,22,0,0), 'sample_t_local': datetime(2023,6,5,13,55,0)}
stations['BANG42'] = {'lat': 47.723367,  'lon': -122.74735, 'sample_time0': '1410-1424', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,22,0,0), 'sample_t_local': datetime(2023,6,5,14,10,0)}
stations['BANG43'] = {'lat': 47.71923987, 'lon': -122.7524348, 'sample_time0': '1430-1447', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,22,30,0), 'sample_t_local': datetime(2023,6,5,14,30,0)}
stations['BANG09'] = {'lat': 47.74300325,  'lon': -122.749475, 'sample_time0': '1500-1511', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,23,0,0), 'sample_t_local': datetime(2023,6,5,15,0,0)}
stations['BANG11'] = {'lat': 47.76017309,  'lon': -122.7339039, 'sample_time0': '1520-1531', 
                      'conc': np.nan, 'sd':np.nan,
                      'model_time': datetime(2023,6,5,23,30,0), 'sample_t_local': datetime(2023,6,5,15,20,0)}
lat = []; lon=[]; all_sites=[]
for key in stations.keys():
    lat.append(stations[key]['lat']);  lon.append(stations[key]['lon'])
    all_sites.append(key)
    
# particle release location
lon_ESP = -122.729553656; lat_ESP = 47.74284803
# distance of station to particle release location
dist_km = np.zeros(len(lat))
for i in range(len(lat)):
    xx, yy = zfun.ll2xy(lon[i], lat[i], lon_ESP, lat_ESP)
    dist_km[i] = np.sqrt(xx**2 + yy**2)/1000
#%% load eDNA concentration from Ryan
tmp = pd.read_csv('concentration_siteMeans.csv')
eDNA_ddPCR = []; eDNA_ddPCR_std = []
# add eDNA conc to stations dictonary
for key in stations.keys():
    stations[key]['conc'] = tmp.loc[tmp['Site']==key, 'SiteMean'].values[0]
    stations[key]['sd']   = tmp.loc[tmp['Site']==key, 'SiteSD'].values[0]
    eDNA_ddPCR.append(stations[key]['conc'])
    eDNA_ddPCR_std.append(stations[key]['sd'])   
eDNA_ddPCR = np.array(eDNA_ddPCR);  eDNA_ddPCR_std = np.array(eDNA_ddPCR_std)
eDNA_ddPCR_ma = np.ma.masked_where(eDNA_ddPCR==0, eDNA_ddPCR) # mask 0 conc
#%% plot map of sampling station 
mpl.rcParams['axes.linewidth'] = 1
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111);
# topography
tfn = 'PS_10m.nc'
tds = xr.open_dataset(tfn)
tx = tds['lon'][::3].values
ty = tds['lat'][::3].values
tz = tds['z'][::3,::3].values
plt.pcolormesh(tx,ty,tz, cmap=cmocean.cm.topo, shading='nearest', vmin=-150, vmax=150, alpha=0.7)
#% bathymetry contour
CS = plt.contour(lonr, latr, hm2, [20,40,60,80,100,120,140], linewidths=0.2, colors='k',linestyles='--',inline_spacings=-8,)
for c in CS.collections:
    c.set_dashes([(0, (10, 10))])
plt.clabel(CS, levels=[20,40,60, 80, 100, 120, 140],inline_spacing=-8, colors='gray', fontsize=7)
# eDNA conc
plt.scatter(lon, lat, s=eDNA_ddPCR*3, facecolors='r', alpha=0.5, linewidths=1.2, edgecolors='r', linestyles='solid')
# eDNA_conc = 0
for ii in range(len(eDNA_ddPCR)):
    if eDNA_ddPCR[ii] == 0:
        plt.plot(lon[ii], lat[ii], 's', markerfacecolor='w', alpha=0.7, ms=10, linestyle='--', markeredgecolor='k',markeredgewidth=1.0)
# add particle release location
plt.plot(lon_ESP, lat_ESP, '*',c='yellow', markeredgecolor='k', markersize=15)
# sampling trajectory   
for i in range(len(lon)-1):
    point1 = (lon[i], lat[i]);  point2 = (lon[i+1], lat[i+1])
    plt.annotate('', xy=point2, xytext=point1, 
                 arrowprops=dict(arrowstyle='->, head_width=0.3', color='k', lw=0.8, ls='-'))
# annotate station ID    
for i in range(len(lon)):
    plt.text(lon[i]+0.001,lat[i],'B'+all_sites[i][-2:], fontsize=11,color='b', rotation=15)
# plot navy ADCP site: # ADCP_IN
plt.plot(-122.7289599, 47.74464362, 'b^')
plt.plot(-122.7289599, 47.74464362, 'r.', markersize=2)
plt.text(-122.7289599, 47.74464362, 'ADCP(2007)', fontsize=14, color='b', rotation=18)
# add a 1km scale bar
plt.hlines(y=47.76, xmin=-122.7135, xmax=-122.7, color='w')
plt.vlines(x=-122.7135, ymin=47.7592, ymax=47.7608, color='w')
plt.vlines(x=-122.7, ymin=47.7592, ymax=47.7608, color='w')
plt.text(-122.71, 47.7605, '1 km', fontsize=14, c='w')
plt.xticks(ticks=[-122.76,-122.74,-122.72,-122.70], fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Longitude', fontsize=15)
plt.ylabel('Latitude', fontsize=15)
aa = [-122.77, -122.695, 47.715, 47.78]
ax.axis(aa); pfun.dar(ax)

# insert zeta
mpl.rcParams['axes.linewidth'] = 0.2
axin = ax.inset_axes([0.52,0.08,0.45,0.12])
axin.patch.set_alpha(0.8)
ds = xr.open_dataset('ESP_2023.06.01_2023.06.10.nc')
ot = ds.ocean_time.values
hr = [np.datetime64(ot1,'h').astype(int)%24 for ot1 in ot]
zeta = ds.zeta.values
axin.plot(ot, zeta, 'k', lw=1); axin.plot(ot, zeta, 'g.')
axin.set_xticks(ot[::2], color='k')
axin.set_xticklabels(hr[::2], fontsize=10, color='k')
axin.set_xlim(datetime(2023,6,5,12,0,0), datetime(2023,6,6,2))
axin.tick_params(axis='x', rotation=0, color='k')
axin.tick_params(axis='y', labelsize=10, colors='k')
axin.set_xlabel('Hours from 2023.06.05 (UTC)', fontsize=10)
axin.text(0.69,0.145,'Surface\nelevation (m)', fontsize=10, transform=ax.transAxes)
axin.set_ylim(-3,2)
# add a box
xx1 = datetime(2023,6,5,18,51,0)
xx2 = datetime(2023,6,5,23,20,0)
yy1 = 2;  yy2 = -3
# Create a Rectangle patch
rect = patches.Rectangle((xx1, yy1), xx2-xx1, yy2-yy1, linewidth=1, edgecolor='r', facecolor='r', alpha=0.1)
axin.add_patch(rect)

# insert wind
axin = ax.inset_axes([0.52,0.235,0.45,0.12])
axin.patch.set_alpha(0.8)
fn = 'zeta_wind_2023.05.24_hc11_v01_uu0k.p'
dict_tmp = pickle.load(open(fn,'rb'))    
zeta = dict_tmp['zeta']
Uwind = dict_tmp['Uwind']
Vwind = dict_tmp['Vwind']
t_tmp = dict_tmp['t']
t = [t1[0] for t1 in t_tmp]
x_pos = np.arange(len(Uwind))
y_pos = np.zeros(len(Vwind))
axin.quiver(t,y_pos,Uwind,Vwind,scale_units='y', scale=1.0,width=0.003, headwidth=2)
axin.set_ylim([-10,10])
axin.plot(t,y_pos, c='k', lw=0.1)
axin.tick_params(axis='y', labelsize=10, colors='k')
hr = [np.datetime64(t1,'h').astype(int)%24 for t1 in t]
axin.set_xticks([datetime(2023,6,3), datetime(2023,6,5), datetime(2023,6,6), datetime(2023,6,8)], color='w')
axin.set_xticklabels(['06/03', '06/05', '06/06', '06/08'], fontsize=9, color='k')
axin.set_xlim(datetime(2023,6,3,0), datetime(2023,6,8,0))
axin.text(0.53,0.33,'Wind (m/s)', fontsize=10, transform=ax.transAxes)
xx1 = datetime(2023,6,5,18,51,0)
xx2 = datetime(2023,6,5,23,20,0)
yy1 = 10; yy2 = -10
# Create a Rectangle patch
rect = patches.Rectangle((xx1, yy1), xx2-xx1, yy2-yy1, linewidth=1, edgecolor='r', facecolor='r', alpha=0.1)
axin.add_patch(rect)
#plt.savefig('June2023_sample_location', dpi=300, bbox_inches='tight')
