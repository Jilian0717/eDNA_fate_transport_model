"""
count particles in each model grid cell, 3D
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import os
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
#from datetime import datetime, timedelta
import datetime as dt
import cmocean.cm as cmo
import pickle
from datetime import datetime, timedelta
import sys

# grid file
grd='grid.nc'
dsg = xr.open_dataset(grd)
latr = dsg['lat_rho'].values
lonr = dsg['lon_rho'].values
h = dsg['h'].values
maskr = dsg['mask_rho'].values
hm = np.ma.masked_where(maskr==1, h)
dx=1/dsg.pm.values  # m
dy=1/dsg.pn.values # m
area = dx*dy

# load delta pier release
rel_dir0 = '../release/'; 
rel_dir = Path(rel_dir0)
fn_list = sorted(rel_dir.glob('release_*'))

AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
      dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
# get grid info
G, S, T = zrfun.get_basic_info('/sciclone/scr-lst/jxiong/LO_roms/hc11_v01_uu0k/f2023.06.05/ocean_his_0020.nc')

continous_cnt = 50 # including previous 50 hours release
target_date = datetime(2023,6,6,3,0,0) # utc

# find initial file
initial_file_time = target_date - timedelta(hours=continous_cnt)
hour0 = initial_file_time.hour
day0 = initial_file_time.day
mon0 = initial_file_time.month
yr0 = initial_file_time.year
initial_file = Path(rel_dir0+'/release_'+str(yr0)+'.'+('00'+str(mon0))[-2:]+'.'+('00'+str(day0))[-2:]+'_sh'+('00'+str(hour0))[-2:]+'.nc')
if initial_file in fn_list:
    ic = fn_list.index(initial_file)

minute = target_date.minute
jc = ic + continous_cnt  # including 50 hourly release

ifile0 = 0
aa,bb = latr.shape; kk=30
counts = np.zeros([aa,bb,kk,len(fn_list)]);  
counts_age = np.zeros([aa,bb,kk,len(fn_list)]);
counts_conc = np.zeros([aa,bb,kk,len(fn_list)]);
    
for ifile in range(ic,jc+1):
    fn = fn_list[ifile]; print(fn)
    sys.stdout.flush()
    dsr = xr.open_dataset(fn, decode_times=False)
    lon_par = dsr.lon.values
    lat_par = dsr.lat.values
    z_par = dsr.z.values
    zeta = dsr.zeta.values
    h_par = dsr.h.values

    # get a list of datetimes
    ot0 = dsr.ot.values
    dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot0]
    # make a mask that is False from the time a particle first leaves the domain
    # and onwards
    ib_mask = np.ones(lon_par.shape, dtype=bool)
    ib_mask[lon_par < AA[0]] = False
    ib_mask[lon_par > AA[1]] = False
    ib_mask[lat_par < AA[2]] = False
    ib_mask[lat_par > AA[3]] = False
    NTS, NPS = lon_par.shape
    for pp in range(NPS):
        tt = np.argwhere(ib_mask[:,pp]==False)
        if len(tt) > 0:
            ib_mask[tt[0][0]:, pp] = False

    # and apply the mask to lon and lat
    lon_par_ma = lon_par.copy()  # masked
    lat_par_ma = lat_par.copy()  
    z_par_ma = z_par.copy()
    lon_par_ma[~ib_mask] = np.nan
    lat_par_ma[~ib_mask] = np.nan
    z_par_ma[~ib_mask] = np.nan

    #% find the nearest grid cell indices for particles
    latr0 = latr[:,0]; lonr0 = lonr[0,:]
    ilat = np.zeros([NPS])
    ilon = np.zeros([NPS]) 
    ik   = np.zeros([NPS])
    
    dt_array = np.array(dt_list)
    # find the closest time
    nti = np.where(np.abs(dt_array-target_date) == np.abs(dt_array-target_date).min())[0][0] 
      
    for npi in range(NPS):
        ilat[npi] = zfun.find_nearest_ind(latr0, lat_par[nti,npi])
        ilon[npi] = zfun.find_nearest_ind(lonr0, lon_par[nti,npi]) 
        z_rho = zrfun.get_z(np.array(h_par[nti, npi]), np.array(zeta[nti, npi]), S, only_rho=True, only_w=False) 
        ik[npi]  = zfun.find_nearest_ind(z_rho, z_par_ma[nti, npi])
              
    ilat = ilat.astype(int); ilon = ilon.astype(int); ik = ik.astype(int)
      
    for npi in range(NPS): 
        if not np.isnan(lon_par_ma[nti, npi]):  # only count particles inside the domain
            counts[ilat[npi], ilon[npi], ik[npi], ifile] += 1
            counts_age[ilat[npi], ilon[npi], ik[npi], ifile] = nti/2 # in hour
            #counts_conc[ilat[npi], ilon[npi], ik[npi], ifile] += C0*np.exp(-decay*nti)
    ifile0 += 1           
#%% 
counts_mean_age = (np.sum(counts*counts_age, axis=3))/(np.sum(counts, axis=3)) # weighted by particle #
counts_mean_age = np.ma.masked_where(counts_mean_age==0, counts_mean_age)

#counts_mean_conc = np.sum(counts_conc, axis=2)
#counts_mean_conc = np.ma.masked_where(counts_mean_conc==0, counts_mean_conc)
counts_all = np.sum(counts, axis=3)
        
#% save plotting variables
tmp = {"counts_all": counts_all,
 #      "counts_mean_conc": counts_mean_conc,
       "counts_mean_age": counts_mean_age}  #
pickle.dump(tmp, open(str(fn)[11:-3].replace('.','_')+"{:02d}".format(minute)+'_3D_previous_'+str(continous_cnt)+'hrs'+".p", "wb"))
