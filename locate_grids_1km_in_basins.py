#!/home/jmframe/programs/anaconda3/bin/python3
# jmframe: I am making this script to loop through the nwm land surface files (LDAS)
# and isolate the values within a shapefile. Save their mean values for LSTM
import datetime
import geopandas as gpd
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
from osgeo import gdal
import pandas as pd
import pickle
import rioxarray
import shapely
from shapely.geometry import Polygon, Point, MultiPolygon
import shapefile
import sys 
import time
import tools
import xarray as xr

# locate the data
dataName = 'camels_all.csv'
#VdataLoc = dataDir + dataName
# load the data with pandas
pd_camatt = pd.read_csv(dataName, sep=',', index_col='gauge_id')

# import shapefile using geopandas
# Comes in with spatial reference 4269
gpd_camels = gpd.read_file('camels_shapes/HCDN_nhru_final_671.shp')
pd_camels = pd.DataFrame(gpd_camels)
pd_camels.set_index('hru_id', inplace=True, drop=True)

# Example file to get the lat/lons for the camels basins
fname = 'LDASOUT/199301010600.LDASOUT_DOMAIN1.comp'
gname = 'geo_em.d01_1km.nc'

# Get the info in numpy array. Especially the lat,lon, which will be transformed
np_ldas = tools.LoadLDAS(fname)
# Example file with data
nc_ldas = nc.Dataset(fname)
# Example file to get the lat/lons for the camels basins
nc_geo_em = nc.Dataset(gname)

# Get the latitude and longitudes
glat = nc_geo_em['XLAT_M'][0,:,:]
glon = nc_geo_em['XLONG_M'][0,:,:]

# Names of the features in this dataset.
featNam = ['ACCET','FIRA','FSA','FSNO','HFX','LH','SNEQV','SNOWH',\
           'SOIL_M1','SOIL_M2','SOIL_M3','SOIL_M4',\
           'SOIL_W1','SOIL_W2','SOIL_W3','SOIL_W4',\
           'TRAD','UGDRNOFF']
# The netcdf has the soil layers stacked in another dimension.
ncNam = ['ACCET','FIRA','FSA','FSNO','HFX','LH','SNEQV','SNOWH',\
         'SOIL_M','SOIL_W','TRAD','UGDRNOFF']

# find nearest cells to camels shapefiles
cam_grid_i = {}
cam_grid_v = {}
if False:
    for cam in list(pd_camels.index.values):
        clat = pd_camatt.loc[cam, 'gauge_lat']
        clon = pd_camatt.loc[cam, 'gauge_lon']
        A = np.array(np.sqrt(np.square((glat-clat))+np.square((glon-clon))))
        imin=np.where(A == np.amin(A))
        minlat = imin[0]
        minlon = imin[1]
        cam_grid_i[cam] = [minlat[0],minlon[0]]
        cam_grid_v[cam] = [glat[minlat,minlon][0],glon[minlat,minlon][0]]
        with open('camels_grid_indices.p', 'wb') as handle:
            pickle.dump(cam_grid_i, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open('camels_grid_values.p', 'wb') as handle:
            pickle.dump(cam_grid_v, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    with open('camels_grid_indices.p', 'rb') as pf:
        cam_grid_i = pickle.load(pf)
    with open('camels_grid_values.p', 'rb') as pf:
        pickle.load(pf)


# Some of the basin shapefiles are multi polygons. 
# But when that is the case there is one polygone for the vast majority of the basin.
# So in this loop I am going to identify, when it is a multi polygon, 
# which has the most points.
# Main polygon Index:
mpindx = {new_list: np.nan for new_list in pd_camels.index.values}
# Loop though basins and find the polygon with the largest number of coordinants
for cam in pd_camels.index.values:
    # Should take when multipolygons come up
    try:
        nshapes=len(gpd_camels.loc[cam].geometry)
        max_n = 0
        for j in range(0,nshapes+1):
            npoints = len(gpd_camels.loc[cam].geometry[j].exterior.coords.xy[0])
            if npoints > max_n:
                max_n = npoints
                mpindx[cam] = j
    # This is a single polygone, so there is no need to find the max number of points.
    except:
        pass

# Make a list of the bounding latitude and longitudes
cambnd ={new_list : [] for new_list in pd_camels.index.values}
for cam in list(pd_camels.index.values):
    cambnd[cam] = pd_camels.loc[cam,'geometry'].bounds

# Here I am going to try to loop through each camels basin, \
# then check all the cells within the bounds.
if False:
    # points in Camels Basins
    ncpindx = {new_list : [] for new_list in pd_camels.index.values}
    for cam in list(pd_camels.index.values):
        print(cam)
        ix = cam_grid_i[cam][1]
        iy = cam_grid_i[cam][0]
        for iup in range(iy,glat.shape[0]):
            llat = glat[iup,ix]
            if llat > cambnd[cam][3]:
                maxilat = iup+1
                break
        for idown in reversed(range(iy)): 
            llat = glat[idown,ix]
            if llat < cambnd[cam][1]:
                minilat = idown-1
                break
        for iright in range(ix,glon.shape[1]):
            llon = glon[iy,iright]
            if llon > cambnd[cam][2]:
                maxilon = iright+1
                break
        for ileft in reversed(range(ix)): 
            llon = glon[iy,ileft]
            if llon < cambnd[cam][0]:
                minilon = ileft-1
                break
        print([maxilat,minilat,maxilon,minilon])
        for y in range(minilat,maxilat):
            for x in range(minilon,maxilon):
                llon = glon[y,x]
                llat = glat[y,x]
                # If the example value is null, then don't do anything and continue loop
                if np_ldas[0,y,x] <= -999:
                    pass
                else:
                    point = Point([llon,llat])
                    if mpindx[cam] >= 0: # For a multi-array. 
                        #Use array with most point
                        if gpd_camels.loc[cam].geometry[mpindx[cam]].contains(point):
                            ncpindx[cam].append([x,y])
                            break
                    else:
                        if gpd_camels.loc[cam].geometry.contains(point):
                            ncpindx[cam].append([x,y])
                            break
            with open('nc_camels_locs_cam1.p', 'wb') as pf:
                pickle.dump(ncpindx, pf, protocol=pickle.HIGHEST_PROTOCOL)
else:
    with open('nc_camels_locs_cam1.p', 'rb') as pf:
        ncpindx = pickle.load(pf)

# Loop through the files and collect the data we need.

if True:

    # initializing time for data aggregation
    loopTime = datetime.datetime(year=1993,month=1,day=1,hour=6)
    timeDelta = datetime.timedelta(hours=3)
    date_list = []

    # Setting datetime for dataframe
    while loopTime < datetime.datetime(year=1993,month=1,day=3,hour=12):
        date_list.append(loopTime)
        loopTime = loopTime + timeDelta
    loopTime = datetime.datetime(year=1993,month=1,day=1,hour=6)

    dynfeaAve = {}
    # Creating pandas data frame
    # Initialize the data frame with all the dynamic features
    for cam in list(pd_camels.index.values):
        dynfeaAve[cam] = pd.DataFrame(index=date_list, columns=featNam)

    # Main loop through time. At each time step gather the data from the files for each basin
    # Then Save the pickle file after each time, so we can re-start whenever.
    for loopTime in date_list:
        # get the file name for this particulare date and time.
        url, fname = tools.construct_LDAS_name(\
                    loopTime.year,loopTime.month,loopTime.day,loopTime.hour)
        print(fname)
        # check that there is a matching file.
        t1 = time.time()
        pe = os.path.exists(fname)
        print('Time to check data',time.time() - t1)

        # If there is a matching file, then open it up
        if pe: 
            np_ldas = tools.LoadLDAS(fname)
            print('Time to load data',time.time() - t1)

            # Getting the values by SLICING THE CAMELS BASINS, the features and the coords
            for cam in pd_camels.index.values:
                grdX = np.array(ncpindx[cam])[:,1]
                grdY = np.array(ncpindx[cam])[:,0]
                M = np_ldas[:,grdX,grdY]
                for ncol, fea in enumerate(featNam):
                    dynfeaAve[cam].loc[loopTime,fea] = np.mean(M[ncol,:])

            print('Time to extract data',time.time() - t1) 
        # If the file doesn't exist, then place NaNs, 
        # that way each day will have a values no matter what.
        else:
            print('file not found')
            for cam in pd_camels.index.values:
                for ncol, fea in enumerate(featNam):
                    dynfeaAve[cam].loc[loopTime,fea] = np.nan
        # Save the pickle file at each time step, so we can start from where we left off.
        with open('dynamic_features_nwm.p', 'wb') as pf:
            pickle.dump(dynfeaAve, pf, protocol=pickle.HIGHEST_PROTOCOL)
else:
    with open('dynamic_features_nwm.p', 'rb') as pf:
        dynfeaAve = pickle.load(pf)
