import numpy as np
import xarray as xr
import netCDF4 as nc
import urllib.request
import pandas as pd
import geopandas as gpd
from rasterio import features
from affine import Affine

def camels_2_nhdp():

	camels = pd.read_table('camels_name.txt',sep=';')
	camels_id = camels.loc[:,'gauge_id']

	convert_map = pd.read_table('usgs_nhdp_id_map.csv',sep=',')
	nhdp_map = convert_map.loc[:,'NHDp']
	usgs_map = convert_map.loc[:,'USGS']
	
	nCamels = camels_id.shape[0]

	nhdp_id = np.full(nCamels,np.nan,dtype=np.int64)
	for c in range(nCamels):
		idex = np.where(camels_id[c] == usgs_map)
		assert(len(idex)==1)
		idex = idex[0]
		nhdp_id[c] = nhdp_map[idex]

	return camels_id, nhdp_id


def camels_nwm_indexes(nhdp_id):

	url = 'https://griffin-objstore.opensciencedatacloud.org/nwm-archive/2014/201409080000.CHRTOUT_DOMAIN1.comp'
	urllib.request.urlretrieve(url, 'temp.nc') 

	dataset = nc.Dataset('temp.nc')
	qqs = dataset.variables['streamflow'][:]
	ids = dataset.variables['feature_id'][:]
	dataset.close()

	idx = np.full(len(nhdp_id),np.nan,dtype=np.int)
	for s in range(len(nhdp_id)):
		itemp = np.where(nhdp_id[s] == ids)
		assert(len(itemp)==1)
		idx[s] = itemp[0]

	return idx


def load_nwm_url(url,fname,site_dex,check_ids):

	nSites = len(site_dex)

	urllib.request.urlretrieve(url,fname) 
	dataset = nc.Dataset(fname)
	qqs = dataset.variables['streamflow'][:]
	ids = dataset.variables['feature_id'][:]
	dataset.close()

	site_qqs = qqs[site_dex]
	site_ids = ids[site_dex]

	assert((site_ids == check_ids).all)

	return site_qqs


def construct_CHRT_name(y,m,d,h):

	ystr = str(y)
	mstr = str(m)
	dstr = str(d)
	hstr = str(h)

	if m < 10: 
		mstr = '0'+mstr

	if d < 10: 
		dstr = '0'+dstr

	if h < 10: 
		hstr = '0'+hstr

	fname = ystr+mstr+dstr+hstr+'00.CHRTOUT_DOMAIN1.comp'
	url = 'https://griffin-objstore.opensciencedatacloud.org/nwm-archive/'+ystr+'/'+fname
	fname = 'CHRTOUT/'+fname

	return url, fname


def construct_LDAS_name(y,m,d,h):

	ystr = str(y)
	mstr = str(m)
	dstr = str(d)
	hstr = str(h)

	if m < 10: 
		mstr = '0'+mstr

	if d < 10: 
		dstr = '0'+dstr

	if h < 10: 
		hstr = '0'+hstr

	fname = ystr+mstr+dstr+hstr+'00.LDASOUT_DOMAIN1.comp'
	url = 'https://griffin-objstore.opensciencedatacloud.org/nwm-archive/'+ystr+'/'+fname
	fname = 'LDASOUT/'+fname

	return url, fname

def construct_LDAS_name_v2(y,m,d,h):
	ystr = str(y)
	mstr = str(m)
	dstr = str(d)
	hstr = str(h)
	if m < 10: 
		mstr = '0'+mstr

	if d < 10: 
		dstr = '0'+dstr

	if h < 10: 
		hstr = '0'+hstr
	fname = ystr+mstr+dstr+hstr+'00.LDASOUT_DOMAIN1.comp'
	fname = '/home/NearingLab/data/nwm/v2/LDAS/'+ystr+'/'+fname
	return fname

def construct_RT_name_v2(y,m,d,h):
	ystr = str(y)
	mstr = str(m)
	dstr = str(d)
	hstr = str(h)
	if m < 10: 
		mstr = '0'+mstr

	if d < 10: 
		dstr = '0'+dstr

	if h < 10: 
		hstr = '0'+hstr
	fname = ystr+mstr+dstr+hstr+'00.RTOUT_DOMAIN1.comp'
	fname = '/home/NearingLab/data/nwm/v2/RT/'+fname
	return fname

def LoadLDAS(fname):
    ldas = nc.Dataset(fname)
    # ACCET Accumulated evapotranspiration (mm)
    ACCET = ldas.variables['ACCET'][:]
    #FIRA: Total net long wave (LW) radiation to atmosphere (W m-2)
    FIRA = ldas.variables['FIRA'][:]
    #FSA: Total absorbed Short Wave (SW) radiation (W m-2)
    FSA = ldas.variables['FSA'][:]
    #FSNO: Snow cover fraction on the ground (fraction)
    FSNO = ldas.variables['FSNO'][:]
    #HFX: Total sensible heat to the atmosphere (W m-2)
    HFX = ldas.variables['HFX'][:]
    #LH: Latent heat to the atmosphere (W m-2)
    LH = ldas.variables['LH'][:]
    #SNEQV: Snow water equivalent (kg m-2)
    SNEQV = ldas.variables['SNEQV'][:]
    #SNOWH: Snow depth (m)
    SNOWH = ldas.variables['SNOWH'][:]
    #SOIL_M: Volumetric soil moisture (m3 m-3)
    SOIL_M0 = ldas.variables['SOIL_M'][:,:,0,:]
    #SOIL_M: Volumetric soil moisture (m3 m-3)
    SOIL_M1 = ldas.variables['SOIL_M'][:,:,1,:]
    #SOIL_M: Volumetric soil moisture (m3 m-3)
    SOIL_M2 = ldas.variables['SOIL_M'][:,:,2,:]
    #SOIL_M: Volumetric soil moisture (m3 m-3)
    SOIL_M3 = ldas.variables['SOIL_M'][:,:,3,:]
    #SOIL_W: Liquid volumetric soil moisture (m3 m-3)
    SOIL_W0 = ldas.variables['SOIL_W'][:,:,0,:]
    #SOIL_W: Liquid volumetric soil moisture (m3 m-3)
    SOIL_W1 = ldas.variables['SOIL_W'][:,:,1,:]
    #SOIL_W: Liquid volumetric soil moisture (m3 m-3)
    SOIL_W2 = ldas.variables['SOIL_W'][:,:,2,:]
    #SOIL_W: Liquid volumetric soil moisture (m3 m-3)
    SOIL_W3 = ldas.variables['SOIL_W'][:,:,3,:]
    #TRAD: Surface radiative temperature (K)
    TRAD = ldas.variables['TRAD'][:]
    #UGDRNOFF: Accumulated underground runoff (mm)
    UGDRNOFF = ldas.variables['UGDRNOFF'][:]
    np_ldas = np.array([ACCET,FIRA,FSA,FSNO,HFX,LH,SNEQV,SNOWH,\
              SOIL_M0,SOIL_M1,SOIL_M2,SOIL_M3,\
              SOIL_W0,SOIL_W1,SOIL_W2,SOIL_W3,\
              TRAD,UGDRNOFF])
    np_ldas = np_ldas[:,0,:,:]
    return np_ldas 

def LoadRT(fname):
    rt = nc.Dataset(fname)
    sfcheadsubrt = rt.variables['sfcheadsubrt'][:]
    zwattablrt = rt.variables['zwattablrt'][:]
    np_rt = np.array([sfcheadsubrt, zwattablrt])
    return np_rt

