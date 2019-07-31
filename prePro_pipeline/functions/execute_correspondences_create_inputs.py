#2D pipeline
import numpy as np
from skimage.morphology import label
from scipy.spatial import cKDTree as KDTree
import pandas as pd
import itertools
from tqdm import tqdm

def Execute_Correspondences_CreateInputs(candidates,normalized_images,im_th,cycle,channels,nbit):
    inputs_df=pd.DataFrame(columns=['cycle','ch','x','y','Intensities_window_5x5'])
    max_df=pd.DataFrame(columns=['I_T','I_G','I_C','I_A','x_T','y_T','x_G','y_G','x_C','y_C','x_A','y_A','cycle'])
    
    cc, n_c = label(np.amax(candidates[cycle,2:channels,:,:],axis=0),return_num=True,connectivity=1)
    conn_components = np.zeros((4,candidates.shape[-2],candidates.shape[-1]))
    for ch in range(4):
        conn_components[ch,:,:] = np.multiply(cc,candidates[cycle,ch+2,:,:])
        
    for i in tqdm(range(1,n_c+1)):
        ch,y,x = np.where(conn_components==i)
        kdT_tmp = KDTree(np.array([x,y]).T)
        if len(list(itertools.combinations(np.arange(len(x)),2)))==len(kdT_tmp.query_pairs(2,p=1)):   # if connected components is too large (likely cover more signals) then split it
            df=pd.Series(data={ 'I_T':np.nan,'I_G':np.nan,'I_C':np.nan,'I_A':np.nan,'x_T':np.nan,'y_T':np.nan,'x_G':np.nan,'y_G':np.nan,'x_C':np.nan,'y_C':np.nan,'x_A':np.nan,'y_A':np.nan,'cycle':cycle})
            df=df[['I_T','I_G','I_C','I_A','x_T','y_T','x_G','y_G','x_C','y_C','x_A','y_A','cycle']]  
            for j in range(len(x)):         
                df.iloc[ch[j]] = im_th[cycle,ch[j]+2,y[j],x[j]]
                df.iloc[ch[j]*2+4]= x[j]
                df.iloc[ch[j]*2+4+1]= y[j]
            I=df['I_T':'I_A']
            col=I[I==np.nanmax(I)].index[0] #retrieving the column 
            tomove=df.index.get_loc(col) #column index to reach the correct columns coordinates
            x_ch=int(df[tomove*2+4])
            y_ch=int(df[tomove*2+4+1])
            ch_idx=tomove
            cycle=int(df['cycle'])
            rect=normalized_images[cycle,ch_idx+2,y_ch-2:y_ch+3,x_ch-2:x_ch+3] 
            if not rect.size==0:
                rect=(rect-np.amin(rect))/(np.amax(rect)-np.amin(rect))
                rect=rect-np.mean(rect)
            row=pd.Series(data={'cycle':cycle,'ch':ch_idx+2,'x':x_ch,'y':y_ch,'Intensities_window_5x5':rect})
                        
            inputs_df=inputs_df.append(row,ignore_index=True)
            max_df=max_df.append(df,ignore_index=True)
        else:
            coords = np.vstack((x,y))
            coords_unique = np.unique(coords,axis=1)
            for j in range(coords_unique.shape[-1]):
                coords_tmp = coords_unique[:,j][:, np.newaxis]
                coords_idx = np.argwhere(np.all(coords==coords_tmp,axis=0)).reshape((-1,))
                df=pd.Series(data={ 'I_T':np.nan,'I_G':np.nan,'I_C':np.nan,'I_A':np.nan,'x_T':np.nan,'y_T':np.nan,'x_G':np.nan,'y_G':np.nan,'x_C':np.nan,'y_C':np.nan,'x_A':np.nan,'y_A':np.nan,'cycle':cycle})
                df=df[['I_T','I_G','I_C','I_A','x_T','y_T','x_G','y_G','x_C','y_C','x_A','y_A','cycle']] 
                for k in coords_idx:
                    df.iloc[ch[k]] = im_th[cycle,ch[k]+2,y[k],x[k]]
                    df.iloc[ch[k]*2+4]= x[k]
                    df.iloc[ch[k]*2+4+1]= y[k]
                I=df['I_T':'I_A']
                col=I[I==np.nanmax(I)].index[0] #retrieving the column 
                tomove=df.index.get_loc(col) #column index to reach the correct columns coordinates
                x_ch=int(df[tomove*2+4])
                y_ch=int(df[tomove*2+4+1])
                ch_idx=tomove
                cycle=int(df['cycle'])
                rect=normalized_images[cycle,ch_idx+2,y_ch-2:y_ch+3,x_ch-2:x_ch+3] 
                if not rect.size==0:
                    rect=(rect-np.amin(rect))/(np.amax(rect)-np.amin(rect))
                    rect=rect-np.mean(rect)
                row=pd.Series(data={'cycle':cycle,'ch':ch_idx+2,'x':x_ch,'y':y_ch,'Intensities_window_5x5':rect})

                inputs_df=inputs_df.append(row,ignore_index=True)
                max_df=max_df.append(df,ignore_index=True)

            
    return {'max_df':max_df, 'inputs_df':inputs_df}
