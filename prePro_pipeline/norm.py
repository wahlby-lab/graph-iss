""" Compute normalization intervals for each channel and round.
    
    Parameters
    ----------
    sys.argv[1] : input csv file contaning an array of input images
    sys.argv[2] : upper percentile value of 98th percentile distri-
        bution of image patches for signal level estimation.
    sys.argv[3] : number of running threads. Each thread run over 
        an image.
    sys.argv[4] : output csv file where to store computed intervals
    sys.argv[5] : number of random patches used to estimate norma-
        lization intervals
"""

import numpy as np
from skimage import io
from scipy import stats
import sys
import pandas as pd
import pickle
from joblib import Parallel,delayed
from sklearn.feature_extraction.image import extract_patches_2d 

def faster_mode1D(a):
    arr = np.asarray(a) # would be _chk_array
    v, c = stats.find_repeats(arr)
    if len(c) == 0:
        arr.sort() # mimic first value behavior
        return arr[0], 1.
    else:
        pos = c.argmax()
        return v[pos], c[pos]

def runParallel(row,seed):
    out = []
    img = io.imread(row.File)
    patch_size = 128 
    patches = extract_patches_2d(img, (patch_size,patch_size), max_patches=int(sys.argv[5]), random_state=int(seed))
    del(img)
    patch = []
    nonZero_px = patch_size*patch_size/4*3
    for i in range(len(patches)):
        if len(patches[i][patches[i]!=0]>=nonZero_px): #(128 x 128 x 3)
            patch.append(patches[i])
    del(patches)

    # Lower bound
    bkg = np.mean([faster_mode1D(patch[i][patch[i]!=0])[0] for i in range(len(patch))])
    print(bkg)
    out.append(bkg)

    # Upper bound
    signal = np.percentile([np.percentile(patch[i],float(sys.argv[2])) for i in range(len(patch))],98)
    out.append(signal)

    return out


imgCSV = pd.read_csv(sys.argv[1],sep='\t')
n_chs = 6
n_cycles = int(len(imgCSV)/6)
seed = np.random.randint(1,100)
res = Parallel(n_jobs=int(sys.argv[3]))(delayed(runParallel)(row,seed) for i, row in imgCSV.iterrows())

img_stats = np.zeros((n_cycles, n_chs, 2))
for i in range(0,len(res),6):
    img_stats[int(i/n_chs),5,0] = res[i+5][0]; img_stats[int(i/n_chs),5,1] = res[i+5][1] # chanA
    img_stats[int(i/n_chs),4,0] = res[i+4][0]; img_stats[int(i/n_chs),4,1] = res[i+4][1]	 # chanC
    img_stats[int(i/n_chs),0,0] = res[i][0]; img_stats[int(i/n_chs),0,1] = res[i][1]	 # chanDO
    img_stats[int(i/n_chs),3,0] = res[i+3][0]; img_stats[int(i/n_chs),3,1] = res[i+3][1]	 # chanG
    img_stats[int(i/n_chs),1,0] = res[i+1][0]; img_stats[int(i/n_chs),1,1] = res[i+1][1]	 # chanNuclei
    img_stats[int(i/n_chs),2,0] = res[i+2][0]; img_stats[int(i/n_chs),2,1] = res[i+2][1]	 # chanT

pickle.dump(img_stats,open(sys.argv[4],'wb'))
