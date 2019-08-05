import pandas as pd
import numpy as np
import networkx as nx
from .pgm.max_flow_min_cost import runMaxFlowMinCost
from scipy.ndimage.filters import gaussian_filter
from scipy.spatial import cKDTree as KDTree

def signalDecoding(data, d_th, k1, tagList, n_threads, dth_max, prior):   
    print("Starting signal decoding......") 
    data = data.rename(index=str, columns={"cycle":"hyb", "Intensities_window_5x5":"Image", "prob_DNN":"p1", "x":"x", "y":"y", "z":"z"})

    data["p0"] = 1-data["p1"]
    data["hyb"] = data["hyb"] + 1
    Imax = []
    Imax_gf = []
    for i,row in data.iterrows():
        if np.isnan(data.Image[i]).any()==True:
            data = data.drop([i])
        else:
            Imax.append(np.amax(data.Image[i]))
            # Gaussian Smoothing Max
            tmp_gfilter = gaussian_filter(data.Image[i], sigma=0.5, truncate=1.5)
            Imax_gf.append(tmp_gfilter[2,2])
    data['Imax'] = Imax
    data['Imax_gf'] = Imax_gf

    data = data.reset_index(drop=True)
    data = data[['hyb','ch','Image','Imax','Imax_gf','p1','p0','x','y','z']]
     
    res = runMaxFlowMinCost(data, d_th, k1, tagList, n_threads, dth_max, prior)
    
    return res


def basecallingMaxFlowMinCost(G, Dvar, Tvar, max_df, final_norm, DNN_model, num_hybs):    
    barcodes = []
    T_Costs=[]
    D_Costs=[]
    signals_df = []
    for c in nx.connected_components(G):
        c=np.array(list(c))
        c=c[c<=Dvar.X_idx.max()]
        if len(Dvar[(Dvar.X_idx.isin(c))])==num_hybs:
            barcodes.append(np.array(np.hstack((np.matrix(Dvar[(Dvar.X_idx.isin(c))].iloc[0,:][['x','y','z']].values),np.matrix(Dvar[(Dvar.X_idx.isin(c))]['ch'])))).flatten().tolist())
            s_df = pd.DataFrame(data=Dvar[(Dvar.X_idx.isin(c))][['ch','hyb','x','y','z']])
            signals_df.append(s_df)

            k1 = KDTree(s_df[['x','y','z']].values)
            T_Costs.append(np.amax(list(k1.sparse_distance_matrix(k1,np.inf).values())))

            prob_s=[]
            I_s=[]
            for s, row in s_df.iterrows():
                ch = row.ch-2
                if ch==0:
                    #T
                    max_df_tmp = max_df[(max_df.x_T==row.x) & (max_df.y_T==row.y) & (max_df.z_T==row.z) & (max_df.cycle==row.hyb-1)]
                elif ch==1:
                    #G
                    max_df_tmp = max_df[(max_df.x_G==row.x) & (max_df.y_G==row.y) & (max_df.z_G==row.z) & (max_df.cycle==row.hyb-1)]
                elif ch==2:
                    #C
                    max_df_tmp = max_df[(max_df.x_C==row.x) & (max_df.y_C==row.y) & (max_df.z_C==row.z) & (max_df.cycle==row.hyb-1)]
                elif ch==3:
                    #A
                    max_df_tmp = max_df[(max_df.x_A==row.x) & (max_df.y_A==row.y) & (max_df.z_A==row.z) & (max_df.cycle==row.hyb-1)]
                    
                I = np.array([max_df_tmp.I_T.values[0],max_df_tmp.I_G.values[0],max_df_tmp.I_C.values[0],max_df_tmp.I_A.values[0]])
                noNaN = np.argwhere(~np.isnan(I)).flatten()
                prob_r = []
                for n in range(4):
                    if n in noNaN:
                        x_ch = max_df_tmp.iloc[:,n*3+4].values[0]
                        y_ch = max_df_tmp.iloc[:,n*3+5].values[0]
                        z_ch = max_df_tmp.iloc[:,n*3+6].values[0]
                    else:
                        x_ch = row.x
                        y_ch = row.y
                        z_ch = row.z
                        I[n] = final_norm[int(max_df_tmp.cycle),int(n+2),int(z_ch),int(y_ch),int(x_ch)]
                      
                    rect=final_norm[int(max_df_tmp.cycle),int(n+2),int(z_ch),int(y_ch-2):int(y_ch+3),int(x_ch-2):int(x_ch+3)]
                    if rect.size==25:
                        rect=(rect-np.amin(rect))/(np.amax(rect)-np.amin(rect))
                        rect=rect-np.mean(rect)
                        X_data=rect.reshape(1,5,5,1)
                        #getting probabilities from the DNN
                        prob_r.append(DNN_model.predict_proba(X_data)[:,1][0])
                    else:
                        prob_r.append(0)
                prob_r = np.array(prob_r).reshape(-1)
                prob_r_max = prob_r[int(row.ch)-2]
                I_max = I[int(row.ch)-2]

                I_ch = np.delete(I,int(row.ch)-2)
                prob_ch = np.delete(prob_r, int(row.ch)-2)
 
                I_ch = I_ch[prob_ch>0.5]
                prob_ch = prob_ch[prob_ch>0.5] 
                if len(prob_ch):
                    prob_s.append((I_max*prob_r_max)/(I_max*prob_r_max + np.amax(I_ch*prob_ch)))
                else:
                    prob_s.append(prob_r_max)
                I_s.append(I_max)    
            prob_s = np.array(prob_s)
            I_s = np.array(I_s)
            D_Costs.append(np.sum(prob_s))
        
    barcodes = pd.DataFrame(barcodes)            
    
    for i in range(3,len(barcodes.columns)):
        barcodes.iloc[:,i].replace([2,3,4,5], ['T', 'G', 'C', 'A'], inplace=True)
    
    return {'barcodes':barcodes, 'T_Costs': T_Costs, 'D_Costs': D_Costs, 'signals_df': signals_df}

