#!/usr/bin/env python

import xarray as xr
import pandas as pd
import numpy as np
import pickle as pkl
import sys
from pathlib import Path

#start_index, stop_index, year= int(sys.argv[1]), int(sys.argv[2]), str(sys.argv[3])
year = str(sys.argv[1])
arbitrary_start = int(sys.argv[2]) if len(sys.argv) > 2 else None
arbitrary_stop = int(sys.argv[3]) if len(sys.argv) > 3 else None

#HOME_DIR = Path('/home/aaraney/gh/nwm_lstm/data')
HOME_DIR = Path('/home/NearingLab/data/nwm/v2/')
DATA_DIR = Path('/home/NearingLab/data/nwm/v2/CHRT/')

# List of files names
FILES = list((DATA_DIR/ year).glob('*'))
FILES.sort()

# File that contains all of the file names
# with open(HOME_DIR / 'chrt_out_pkl.p', 'rb') as f:
    # domain_files = pkl.load(f)

# Account for going past all indexes
#if stop_index > len(domain_files):
#    stop_index = len(domain_files)

# Exclusive stop
#domain_files = FILES[start_index:stop_index]

camel_site_no = HOME_DIR / 'nwm_camels_basin_id_and_nhd_link.csv'
# camel_site_no = '../data/nearest_com_id_to_missing_camels_simplified.csv'
camel_site_data = pd.read_csv(camel_site_no, dtype={'site_no': str, 'feature_id': int})

camel_site_no = camel_site_data['site_no'].tolist()
camel_site_comid = camel_site_data['feature_id'].tolist()

camel_site_data.set_index('site_no', inplace=True)

# Create a dictionary with site_no's from camel's dataset as keys with empty values
# dynamic_features = dict.fromkeys(list(camel_site_no.site_no))
d_f_keys = pd.DataFrame(columns = ['streamflow', 'q_lateral', 'velocity', 'qSfcLatRunoff', 'qBucket', 'qBtmVertRunoff'])

# Map empty dataframes with desired column headers to each key in dictionary
dynamic_features = {k: d_f_keys for k in camel_site_no}

# df_feature_id = xr.open_dataset(FILES[0]).sel(feature_id = camel_site_comid)

def drop_coords(ds):
    ds = ds.sel(feature_id = camel_site_comid)
    ds = ds.drop_vars(['latitude', 'longitude', 'order', 'crs', 'elevation',])
    ds = ds.drop_dims(['reference_time',]) #'feature_id',])
    return ds.reset_coords(drop=True)


#i = 0
#for _file in range(50, len(domain_files) + 50, 50):

    # Account for the step not including all of the values
    #if _file > len(domain_files):
        #_file = len(domain_files)
if arbitrary_start:
    start = arbitrary_start
    loop_start = start + 50

else:
    start = 0
    loop_start = 50

stop = len(FILES) + 51

for i in range(loop_start,stop,50):
    try:
        if i < len(FILES):
            with xr.open_mfdataset(FILES[start:i],
        		   preprocess=drop_coords,
        		   combine='nested',
        		   concat_dim='time',
        		   decode_cf=False) as df:
        
                df = xr.decode_cf(df)
                #df.coords['feature_id'] = df_feature_id.coords['feature_id']
                
                df = df.to_dataframe()
        else:
            with xr.open_mfdataset(FILES[start:],
        		   preprocess=drop_coords,
        		   combine='nested',
        		   concat_dim='time',
        		   decode_cf=False) as df:
    
                df = xr.decode_cf(df)
                #df.coords['feature_id'] = df_feature_id.coords['feature_id']
                
                df = df.to_dataframe()

        dynamic_features = {k: v.append(df.loc[camel_site_data.loc[k].values[0]]) for (k, v) in dynamic_features.items()}

    except:
        # Check if dict of dfs has empty df's
        if not dynamic_features[list(dynamic_features.keys())[0]].empty:
            with open(HOME_DIR / 'CHRT_OUT' / 'dynamic_features_ic_2_{}.p'.format(year), 'wb') as pf: 
                pkl.dump(dynamic_features, pf)

        with open(HOME_DIR / 'logs' / 'dynamic_features_log_{}.log'.format(year), 'a+') as f:
            for x in FILES[start:i]:
                try:
                    _df = xr.open_dataset(x)
                except:
                   f.write("Couldn't load: {}\n".format(x))
                   f.write("start, stop: {}, {}\n".format(start, i))
                   FILES.pop(x)
        
    # df = df.drop(columns = ['latitude','longitude']) 
    #i = _file + 1
    
    # df = (df[df.index.get_level_values(0)
    #       .isin(camel_site_no['feature_id'])])
    
    # df = (df.join(camel_site_no
    #    .set_index('feature_id'), on='feature_id', how='left')
    #    .set_index('site_no', append=True)
    #    .reorder_levels(['feature_id', 'site_no', 'reference_time', 'time']))
    

    start += 50

# Fix Multi-indexing
#  dynamic_features = {k: v.reindex(pd.MultiIndex.from_tuples(v.index)) for (k, v) in dynamic_features.items()}

# Remove NA dates
# dynamic_features = {k: v[['streamflow', 'q_lateral', 'velocity', 'qSfcLatRunoff', 'qBucket', 'qBtmVertRunoff']].dropna() for (k, v) in dynamic_features.items()}

# Save dynamic_features dictionary of dataframes to a pickle file
with open(HOME_DIR / 'CHRT_OUT' / 'dynamic_features_{}.p'.format(year), 'wb') as pf:
    pkl.dump(dynamic_features, pf)
