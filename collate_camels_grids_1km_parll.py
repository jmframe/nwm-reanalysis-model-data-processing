#!/home/jmframe/programs/anaconda3/bin/python3
import pickle
import pandas as pd
import numpy as np
import os
import time
from joblib import Parallel, delayed
import multiprocessing
# custom
import tools
import sys

#Start year
s_y = sys.argv[1]
pd_s_y = str(int(s_y))+'-01-01 00:00'
#End year
pd_e_y = str(int(s_y)+1)+'-01-01 00:00'

# Load Gauge IDs
basin_file = "basin_list.txt"
with open(basin_file, 'r') as fp:
    basins = fp.readlines()
    basins = [basin.strip() for basin in basins]
    
nBasins = len(basins)

# Names of the features in this dataset.
features = ['ACCET','FIRA','FSA','FSNO','HFX','LH','SNEQV','SNOWH',\
           'SOIL_M1','SOIL_M2','SOIL_M3','SOIL_M4',\
           'SOIL_W1','SOIL_W2','SOIL_W3','SOIL_W4',\
           'TRAD','UGDRNOFF']

nFeatures = len(features)

# initialize time index
dates = pd.date_range(pd_s_y, pd_e_y, freq="180min")

nTimes = len(dates)

with open('nc_camels_locs_cam1.p', 'rb') as pf:
    ncpindx = pickle.load(pf)

def collect_single_timestep(date):#, basins, ncpindx, features):

    # get file name
    fname = tools.construct_LDAS_name_v2(date.year,
                                         date.month,
                                         date.day,
                                         date.hour)
    print(fname)

    # init storage
    timestep_data_np = np.full([len(basins), len(features)], np.nan)    
    
    # check that there is a matching file.
    pe = os.path.exists(fname)
    
    # If there is a matching file, then open it up
    if pe: 

        # benchmarking
        t1 = time.time()
        
        # load netcdf file
        try:
            np_ldas = tools.LoadLDAS(fname)
            t2 = time.time()
            print('Time to load data', t2-t1) 
            
            # Getting the values by SLICING THE CAMELS BASINS, the features and the coords
            for ib, b in enumerate(basins):
                grdX = np.array(ncpindx[int(b)])[:,1]
                grdY = np.array(ncpindx[int(b)])[:,0]
                try:
                    timestep_data_np[ib,:] = np.nanmean(np_ldas[:,grdX,grdY], axis=1)
                except: 
                    timestep_data_np[ib,:] = np.nan
            t3 = time.time()
            print('Time to extract data', t3-t2)
        except: 
            print('failed to open file: ', fname)
    else:
        print('failed to open file: ', fname)
    return timestep_data_np

num_cores = multiprocessing.cpu_count()
print(f"There are {num_cores} cores available")

# Run the data collation in parallel
results_parallel = Parallel(n_jobs=num_cores)(delayed(collect_single_timestep)(dates[t]) for t in range(nTimes))

def make_dictionary(dates,features,basins):
    # Initialize space to save data
    camels_dictionary = {}
    # Creating pandas data frame
    # Initialize the data frame with all the dynamic features
    for b in basins:
        camels_dictionary[b] = pd.DataFrame(index=dates, columns=features)
    return camels_dictionary
camels_dictionary = make_dictionary(dates,features,basins)

for ib, b in enumerate(basins):
    camels_dictionary[b].loc[:] = np.array(results_parallel)[:,ib,:]

with open('/home/NearingLab/data/nwm/v2/LDAS/camels_basins/dynamic_features_nwm_RT_'+str(s_y)+'.p', 'wb') as pf:
    pickle.dump(camels_dictionary, pf, protocol=pickle.HIGHEST_PROTOCOL)
