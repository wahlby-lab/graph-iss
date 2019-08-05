import numpy as np
from skimage.measure import label 
from scipy.spatial import distance
from scipy import stats
from skimage.morphology import extrema,white_tophat,disk,diamond
import pandas as pd

def faster_mode1D(a):
    arr = np.asarray(a)
    v, c = stats.find_repeats(arr)
    if len(c) == 0:
        arr.sort() # mimic first value behavior
        return arr[0], 1.
    else:
        pos = c.argmax()
        return v[pos], c[pos]


def ExtractCandidates(im_norm,h,radius,nbit):
    """extract signal candidates applying h_maxima transform 
        INPUTS:
        im_norm=normalised image,
        h=h_maxima threshold,
        radius=structuring element radius,
        nbit= encoding"""
    
    # Normalized top_hat filtering
    se=disk(radius)
    im=white_tophat(im_norm,se)
        
    #filtering local maxima
    h_maxima=extrema.h_maxima(im, h,selem=diamond(1))
    label_h_max=label(h_maxima,neighbors=4)
    labels=pd.DataFrame(data={'labels':np.sort(label_h_max[np.where(label_h_max!=0)])})
    dup=labels.index[labels.duplicated() == True].tolist() #find duplicates labels (=connected components) 
    
    #splitting connected regions to get only one local maxima    
    max_mask=np.zeros(im.shape)
    max_mask[label_h_max!=0]=np.iinfo(nbit).max
            
    for i in range (len(dup)):
        r,c=np.where(label_h_max==labels.loc[dup[i],'labels']) #find coord of points having the same label
        meanpoint_x=np.mean(c)
        meanpoint_y=np.mean(r)
        dist=[distance.euclidean([meanpoint_y,meanpoint_x],[r[j],c[j]]) for j in range(len(r))]
        ind=dist.index(min(dist))
        r,c=np.delete(r,ind),np.delete(c,ind) #delete values at ind position.
        max_mask[r,c]=0 #set to 0 points != medoid coordinates
             
    return max_mask

