import numpy as np
from skimage.measure import label 
from scipy.spatial import distance
from scipy import stats
from skimage.morphology import extrema,white_tophat,ball
import pandas as pd

def faster_mode1D(a):
    arr = np.asarray(a) # would be _chk_array
    v, c = stats.find_repeats(arr)
    if len(c) == 0:
        arr.sort() # mimic first value behavior
        return arr[0], 1.
    else:
       pos = c.argmax()
       return v[pos], c[pos]


def ExtractCandidates(im_norm,h,radius,nbit,sigma):
    """extract signal candidates applying h_maxima transform 
    INPUTS:
        im_norm=normalised image,
        h=h_maxima threshold,
        radius=structuring element radius,
        nbit= encoding"""
    
    #top_hat filtering
    se=ball(radius)
    im=white_tophat(im_norm,se)
    connectivity=np.array([[[0, 0, 0],[0, 1, 0],[0, 0, 0]],[[0, 1, 0],[1, 1, 1],[0, 1, 0]],[[0, 0, 0],[0, 1, 0],[0, 0, 0]]]) #3D corss  

    
    #filtering local maxima
    h_maxima=extrema.h_maxima(im,h,selem=connectivity)
    label_h_max=label(h_maxima, neighbors=4)
    labels=pd.DataFrame(data={'labels':np.sort(label_h_max[np.where(label_h_max!=0)])})
    dup=labels.index[labels.duplicated() == True].tolist() #find duplicates labels (=connected components) 
    
    #splitting connected regions to get only one local maxima    
    max_mask=np.zeros(im.shape)
    max_mask[label_h_max!=0]=np.iinfo(nbit).max
            
    for i in range (len(dup)):
        z,r,c=np.where(label_h_max==labels.loc[dup[i],'labels']) #find coord of points having the same label
        meanpoint_x=np.mean(c)
        meanpoint_y=np.mean(r)
        meanpoint_z=np.mean(z)
        dist=[distance.euclidean([meanpoint_z,meanpoint_y,meanpoint_x],[z[j],r[j],c[j]]) for j in range(len(r))]
        ind=dist.index(min(dist))
        z,r,c=np.delete(z,ind),np.delete(r,ind),np.delete(c,ind) #delete values at ind position.
        max_mask[z,r,c]=0 #set to 0 points != medoid coordinates
    
##    overlay_h=color.label2rgb(label_h_max, im, alpha=1, bg_label=0, bg_color=None, colors=[(0,1,0)])   
##    plt.figure()
##    plt.imshow(overlay_h,interpolation='none')
         
    return max_mask

