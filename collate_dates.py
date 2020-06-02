
import pickle as pkl
from pathlib import Path

DIR = Path('/home/NearingLab/data/nwm/v2/CHRT_OUT/')

DIR = list(DIR.glob('*.p'))
DIR.sort()

[x for x in DIR if x.ends_with('1999.p')]
