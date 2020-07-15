#!/usr/bin/python3

import pandas as pd
import numpy as np
import pickle as pkl
import sys
from pathlib import Path
from os.path import basename


# Provide pickle file of dictionaries of dataframes as first command
# line arg
file_1 = sys.argv[1]

with open(file_1, 'rb') as f:
    df = pkl.load(f)


def daily_average(df):
   return(df.groupby(pd.Grouper(freq='1D')).mean())

   # df.index.rename(['comid', 'station_id', 'ref_date', 'time'], inplace  =  True)

   # return(df.reset_index().groupby(pd.Grouper(key = 'time', freq='1D')).mean())

means = {k: daily_average(v) for (k, v) in df.items()}

# means = df.groupby(pd.Grouper(freq='1D')).mean()

file_1 = basename(Path(file_1))
# file_1.split('.')[0]

with open(f'{file_1}_daily_means.p', 'wb') as f:
    pkl.dump(means, f)
