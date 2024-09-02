import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.patches as patches
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
import cmocean.cm as cmo
import pickle
from datetime import datetime, timedelta
from matplotlib import cm
import cmaps
from matplotlib.backend_bases import MouseButton

# grid file
grd='grid.nc'
dsg = xr.open_dataset(grd)
latr = dsg['lat_rho'].values
lonr = dsg['lon_rho'].values
h = dsg['h'].values
mask = dsg['mask_rho'].values
hm = np.ma.masked_where(mask==1, h)

#%
# Information of June's sampling
stations = dict()
stations['BANG16'] = {'lat': 47.7743587,  'lon': -122.7072554,
                      'sample_time0': '1051-1057',
                      'conc': 27.28608569, 'sd':np.nan,  # copies/L
                      'model_time':     datetime(2023,6,5,19,0,0),
                      'sample_t_local': datetime(2023,6,5,10,51,0)}
stations['BANG17'] = {'lat': 47.77476498,  'lon': -122.7139617,
                      'sample_time0': '1103-1109',
                      'conc': 0.0, 'sd':23.28702888,
                      'model_time':     datetime(2023,6,5,19,0,0),
                      'sample_t_local': datetime(2023,6,5,11,3,0)}
stations['BANG15'] = {'lat': 47.77071505, 'lon': -122.7132016,
                      'sample_time0': '1118-1129',
                      'conc': 11.36204826, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,19,30,0),
                      'sample_t_local': datetime(2023,6,5,11,18,0)}
stations['BANG14'] = {'lat': 47.77053826, 'lon': -122.7181893,
                      'sample_time0': '1138-1144',
                      'conc': 16.60320373, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,19,30,0),
                      'sample_t_local': datetime(2023,6,5,11,38,0)}
stations['BANG35'] = {'lat': 47.76205,  'lon': -122.727267,
                      'sample_time0': '1155-1202',
                      'conc': 6.236040799, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,20,0,0),
                      'sample_t_local': datetime(2023,6,5,11,55,0)}
stations['BANG36'] = {'lat': 47.758383,  'lon': -122.732183,
                      'sample_time0': '1209-1219',
                      'conc': 0.0, 'sd':1.673358021,
                      'model_time':     datetime(2023,6,5,20,0,0),
                      'sample_t_local': datetime(2023,6,5,12,9,0)}
stations['BANG37'] = {'lat': 47.755,  'lon': -122.736083,
                      'sample_time0': '1225-1235',
                      'conc': 6.689918809,'sd':np.nan,
                      'model_time':     datetime(2023,6,5,20,30,0),
                      'sample_t_local': datetime(2023,6,5,12,25,0)}
stations['BANG38'] = {'lat': 47.747767,  'lon': -122.741083,
                      'sample_time0': '1244-1257',
                      'conc': 0.0,'sd':np.nan,
                      'model_time':     datetime(2023,6,5,20,30,0),
                      'sample_t_local': datetime(2023,6,5,12,44,0)}
stations['BANG39'] = {'lat': 47.742333,   'lon': -122.742417,
                      'sample_time0': '1305-1318',
                      'conc': 2.220256128, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,21,0,0),
                      'sample_t_local': datetime(2023,6,5,13,5,0)}
stations['BANG40'] = {'lat': 47.738267,  'lon': -122.74495,
                      'sample_time0': '1323-1338',
                      'conc': 6.767606932, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,21,30,0),
                      'sample_t_local': datetime(2023,6,5,13,23,0)}
stations['BANG41'] = {'lat': 47.7297,     'lon': -122.749283,
                      'sample_time0': '1355-1405',
                      'conc': 9.952124767, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,22,0,0),
                      'sample_t_local': datetime(2023,6,5,13,55,0)}
stations['BANG42'] = {'lat': 47.723367,  'lon': -122.74735,
                      'sample_time0': '1410-1424',
                      'conc': 0.0, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,22,0,0),
                      'sample_t_local': datetime(2023,6,5,14,10,0)}
stations['BANG43'] = {'lat': 47.71923987, 'lon': -122.7524348,
                      'sample_time0': '1430-1447',
                      'conc': 614.4670593, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,22,30,0),
                      'sample_t_local': datetime(2023,6,5,14,30,0)}
stations['BANG09'] = {'lat': 47.74300325,  'lon': -122.749475,
                      'sample_time0': '1500-1511',
                      'conc': 0.928269305, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,23,0,0),
                      'sample_t_local': datetime(2023,6,5,15,0,0)}
stations['BANG11'] = {'lat': 47.76017309,  'lon': -122.7339039,
                      'sample_time0': '1520-1531',
                      'conc': 13.53247011, 'sd':np.nan,
                      'model_time':     datetime(2023,6,5,23,30,0),
                      'sample_t_local': datetime(2023,6,5,15,20,0)}

# find location of detections
detected_sites = []; lon=[]; lat=[];
for key in stations.keys():
    lon.append(stations[key]['lon'])
    lat.append(stations[key]['lat'])
    if stations[key]['conc']>0:
        detected_sites.append(key)

eDNA_ddPCR = []
eDNA_ddPCR_std = []
for key in stations.keys():
    eDNA_ddPCR.append(stations[key]['conc'])
    eDNA_ddPCR_std.append(stations[key]['sd'])
eDNA_ddPCR = np.array(eDNA_ddPCR)
eDNA_ddPCR_std = np.array(eDNA_ddPCR_std)
        
#%% load vol and zeta
tmp = pd.read_pickle('vol_zeta_2023.05.23_hc11_v01_uu0k.p')
vol_avg = np.transpose(tmp['vol']/tmp['cnt'], axes=[1,2,0])
vol_avg[vol_avg==0] = np.nan
zeta_avg = tmp['zeta']/tmp['cnt']
#%% import 3D particle number
fn0='./data/release_2023_06_05_sh1900_3D_previous_50hrs.p'; shading_sites = ['BANG16']
#fn0='./data/release_2023_06_05_sh1930_3D_previous_50hrs.p'; shading_sites = ['BANG14', 'BANG15']
#fn0='./data/release_2023_06_05_sh2000_3D_previous_50hrs.p'; shading_sites = ['BANG35']
#fn0='./data/release_2023_06_05_sh2030_3D_previous_50hrs.p'; shading_sites = ['BANG37']
#fn0='./data/release_2023_06_05_sh2100_3D_previous_50hrs.p'; shading_sites = ['BANG39']
#fn0='./data/release_2023_06_05_sh2130_3D_previous_50hrs.p'; shading_sites = ['BANG40']
#fn0='./data/release_2023_06_05_sh2200_3D_previous_50hrs.p'; shading_sites = ['BANG41']
#fn0='./data/release_2023_06_05_sh2230_3D_previous_50hrs.p'; shading_sites = ['BANG43']
#fn0='./data/release_2023_06_05_sh2300_3D_previous_50hrs.p'; shading_sites = ['BANG09']
#fn0='./data/release_2023_06_05_sh2330_3D_previous_50hrs.p'; shading_sites = ['BANG11']

tmp = pickle.load(open(fn0, "rb")) 
counts_all = tmp['counts_all'] # [lon,lat,vertical] # summed particle # from previous 50hours
counts_mean_age = tmp['counts_mean_age']   
#% cross channel
dict_tmp = pickle.load(open("HC_along_sect3.p",'rb'))
pts = dict_tmp['pts']
lon_tran = []; lat_tran = []
for i in range(len(pts)):
    lon_tran.append(pts[i][0])
    lat_tran.append(pts[i][1])
lon_tran = np.array(lon_tran)
lat_tran = np.array(lat_tran)

xx = lon_tran; yy = lat_tran
for ii in range(len(xx)-1):
    x0 = xx[ii]
    x1 = xx[ii+1]
    y0 = yy[ii]
    y1 = yy[ii+1]
    nn = 50
    if ii == 0:
        x_e = np.linspace(x0, x1, nn)
        y_e = np.linspace(y0,y1, nn)
    else:
        x_e = np.concatenate((x_e, np.linspace(x0, x1, nn)[1:]))
        y_e = np.concatenate((y_e, np.linspace(y0, y1, nn)[1:]))
        
G, S, T = zrfun.get_basic_info('/Users/jilian/Desktop/LiveOcean/LO_roms/hc11_v01_uu0k/f2023.06.05/ocean_his_0020.nc')
h = G['h']
zw = zrfun.get_z(h, zeta_avg, S, only_w=True)
zrho = zrfun.get_z(h, zeta_avg, S, only_rho=True)

ix_2m = (zrho-zrho[-1,:,:])>=-2  # upper 2m
ix_2m = np.transpose(ix_2m, (1,2,0))

import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 0.5

#%--------------------par density in upper 2m -------------------
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
tmp = counts_all/vol_avg
tmp11 = np.ma.masked_where(~ix_2m, tmp)
tmp0 = np.nanmean(tmp11, axis=2)
tmp0[tmp0==0] = np.nan

cs = plt.pcolormesh(lonr,latr, tmp0,  # upper layer
                    cmap = cmo.matter,             
                    shading='auto',#alpha=0.7, 
                    vmin=0, vmax=0.05) # particle number density
cs1 = plt.pcolormesh(lonr, latr, hm, cmap='gray', alpha=0.1)    

# ESP location Jan2023
lon0 = -122.729553656; lat0 = 47.74284803
plt.plot(lon0, lat0, '*',c='yellow', markeredgecolor='k', markersize=15)
plt.plot(lon, lat,'.',c='b', markersize=5)
# detected station
plt.scatter(lon, lat, s=eDNA_ddPCR*3, facecolors='none', edgecolors='r')

# shading detected location
for site in shading_sites:
    plt.scatter(stations[site]['lon'], stations[site]['lat'], s=stations[site]['conc']*3,  
          facecolors='r', alpha=0.5)
    plt.text(stations[site]['lon'], stations[site]['lat'],'B'+site[-2:],\
             fontsize=18, color='b', rotation=15)

aa = [-122.77, -122.695, 47.715, 47.78]
ax.axis(aa); pfun.dar(ax)
plt.xticks(ticks=[-122.76,-122.72],labels = ['-122.76', '-122.72'], fontsize=25)
plt.yticks(ticks=[47.73, 47.75, 47.77], fontsize=25)
tmp=fn0[-37:-20]
plt.text(0.1,0.94,tmp[13:15]+':'+tmp[15:], ha='left', 
         transform=ax.transAxes, fontsize=25, c='k')

# add colorbar
if fn0[-26:-20] == 'sh1900':
    cbar_ax = fig.add_axes([0.45, 0.2, 0.28, 0.02])
    cbar = plt.colorbar(cs, cax=cbar_ax, orientation='horizontal')#, extend='both')
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_title('par#/m$^3$', fontsize=22) 
    
plt.savefig('./figs/'+'upper2m_'+str(fn0).replace('.','_')[-37:-20], dpi=300, bbox_inches='tight')


