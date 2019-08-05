import pandas as pd
from joblib import Parallel,delayed 
import numpy as np
from skimage.morphology import white_tophat,ball
from functions.execute_correspondences_create_inputs import Execute_Correspondences_CreateInputs

def findCorrespondencesMax_createInputs(candidates,normalized_images,cycles,channels,nbit,n_threads,radius):
    """INPUTS:
        candidates=binary max masks
        normalized_images=normalized images,
        cycles=number of sequencing cycles,
        nbit= encoding"""    
    se=ball(radius)
    im_th = np.zeros(normalized_images.shape).astype(np.float64)
    for cycle in range(cycles):
        for ch in range(2,channels):
            im_th[cycle,ch,:,:,:] = white_tophat(normalized_images[cycle,ch,:,:,:],se)

    candidates[candidates==np.iinfo(nbit).max]=1

    res=Parallel(n_jobs=n_threads)(delayed(Execute_Correspondences_CreateInputs)(candidates, normalized_images, im_th, cycle, channels, nbit) for cycle in range(cycles))

    max_df=pd.concat([res[i]['max_df'] for i in range(len(res))],ignore_index=True)
    inputs_df=pd.concat([res[i]['inputs_df'] for i in range(len(res))],ignore_index=True) 
    inputs_df=inputs_df[['cycle','ch','x','y','z','Intensities_window_5x5']]
    
    return {'inputs_df':inputs_df, 'max_df':max_df}
    
   
