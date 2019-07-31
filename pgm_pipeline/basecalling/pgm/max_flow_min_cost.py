import pandas as pd
import numpy as np
from scipy.spatial import cKDTree as KDTree
import datetime
import networkx as nx
from tqdm import tqdm
import itertools

# Function Gibbs Distribution: convert probablities to energies
def prob2Eng(p):
	return -1.0 * np.log(np.clip(p,0.00001,0.99999))

def runComponent(data, l, d_th, k1, k2, tagList, num_hyb, dth_max, prior):
    if len(data[data.conn_comp ==l]):
        if len(np.unique(data[data.conn_comp ==l].hyb))==len(num_hyb):
            data_tmp = data[data.conn_comp ==l].sort_values(['hyb']).copy()
            data_tmp.reset_index(inplace=True, drop=True)
            Dvar_tmp = pd.DataFrame(data={'x':[],'y':[],'z':[],'hyb':[], 'ch':[], 'X_idx':[]})
            for h in num_hyb:
                for i, row in data_tmp[data_tmp.hyb==h].iterrows():
                    # add detection variables
                    Dvar_tmp = Dvar_tmp.append(pd.Series([row.x,row.y,row.z,h,row.ch,i, prob2Eng(row.p0), prob2Eng(row.p1)],index=['x','y','z','hyb','ch','X_idx','E_0','E_1']),ignore_index=True)
                
            X_idx_tmp=len(Dvar_tmp)
            Tvar_tmp = pd.DataFrame(data={'x_idx':[],'anchestor_x_idx':[], 'descendant_x_idx':[],'E_0':[], 'E_1':[]})
            for h1 in num_hyb[:-1]:
                h2=h1+1
                if len(data_tmp[data_tmp.hyb==h1])==0:
                    continue
                KDTree_h1 = KDTree(data_tmp[data_tmp.hyb==h1][['x','y','z']])
                if len(data_tmp[data_tmp.hyb==h2])==0:
                    continue
        
                KDTree_h2 = KDTree(data_tmp[data_tmp.hyb==h2][['x','y','z']])
                query = KDTree_h1.query_ball_tree(KDTree_h2,dth_max, p=2)
                for i in range(len(query)):
                    if len(query[i]):
                        data_i_idx = data_tmp.index[data_tmp.hyb==h1][i]
                        row_i = data_tmp.loc[[data_i_idx]]
                        for j in range(len(query[i])):
                            data_j_idx = data_tmp.index[data_tmp.hyb==h2][query[i][j]]
                            row_j = data_tmp.loc[[data_j_idx]]
                            d = np.linalg.norm(np.array(row_i[['x','y','z']])-np.array(row_j[['x','y','z']]))
                            # create transition var only if two signals are close enough
                            fI = np.absolute(np.float64(row_i.Imax_gf)-np.float64(row_j.Imax_gf))
                            mu_d = 1/ (1 + k1 * d)
                            mu_I = mu_d / (1 + k2 * fI)
                            Tvar_tmp = Tvar_tmp.append(pd.Series([X_idx_tmp,data_i_idx,data_j_idx, prob2Eng(1-mu_d), prob2Eng(mu_d), mu_d, mu_I],index=['x_idx','anchestor_x_idx', 'descendant_x_idx','E_0', 'E_1', 'mu_D', 'mu_I']),ignore_index=True)
                            X_idx_tmp=X_idx_tmp+1;
        
            if prior=="prior":
                # Exclude paths not in taglist        
                Tvar_tmp.x_idx=np.ravel_multi_index([Tvar_tmp.anchestor_x_idx.astype(np.int),Tvar_tmp.descendant_x_idx.astype(np.int)],(len(Dvar_tmp),len(Dvar_tmp)))
                hyb_sets = [];
                for h in num_hyb:
                    hyb_sets.append(np.unique(np.array(Dvar_tmp[Dvar_tmp.hyb==h].ch)))
                allPaths=list(itertools.product(*hyb_sets))
                allowedPaths_D_list = []
                allowedPaths_T_list = []
                # Remove not existing paths
                if len(allPaths):
                    for tag in range(len(tagList)):
                        if (allPaths == np.array(tagList.iloc[tag,:]).astype(int)).all(1).any():
                            allowedDvar = pd.DataFrame(data={'x':[],'y':[],'z':[],'hyb':[], 'ch':[], 'X_idx':[],'E_0':[],'E_1':[]});
                            for h in num_hyb:
                                allowedDvar = allowedDvar.append(Dvar_tmp[(Dvar_tmp.hyb == h) & (Dvar_tmp.ch==int(tagList.iloc[tag,int(h-1)]))], ignore_index=True)
                            DGraph = nx.Graph()
                            DGraph.add_nodes_from(allowedDvar.X_idx)
                            DGraph.add_edges_from(Tvar_tmp[(Tvar_tmp.anchestor_x_idx.isin(allowedDvar.X_idx)) & (Tvar_tmp.descendant_x_idx.isin(allowedDvar.X_idx))][['anchestor_x_idx', 'descendant_x_idx']].values.tolist())                
                            for c in nx.connected_components(DGraph):
                                if len(allowedDvar[allowedDvar.X_idx.isin(c)].hyb.unique()) == len(num_hyb):
                                    allowedPaths_D_list.append(np.array(list(c)))                      
                                    tmp_dvar = np.array(Tvar_tmp[(Tvar_tmp.anchestor_x_idx.isin(c)) & (Tvar_tmp.descendant_x_idx.isin(c))][['anchestor_x_idx', 'descendant_x_idx']]);
                                    tmp_path = []
                                    for v in range(len(tmp_dvar)):
                                        tmp_path.append(np.ravel_multi_index([int(tmp_dvar[v,0]),int(tmp_dvar[v,1])],(len(Dvar_tmp),len(Dvar_tmp))))
                                    allowedPaths_T_list.append(np.array(tmp_path))
                
                Dvar_tmp = Dvar_tmp[Dvar_tmp.X_idx.isin(np.concatenate(allowedPaths_D_list))]
                Tvar_tmp = Tvar_tmp[Tvar_tmp.x_idx.isin(np.concatenate(allowedPaths_T_list))]
                
                d=dict(zip(Dvar_tmp.X_idx, np.arange(len(Dvar_tmp))))
                Dvar_tmp.X_idx = [d[x] for x in Dvar_tmp.X_idx]
                Tvar_tmp.anchestor_x_idx = [d[x] for x in Tvar_tmp.anchestor_x_idx]
                Tvar_tmp.descendant_x_idx = [d[x] for x in Tvar_tmp.descendant_x_idx]
                
            Dvar_tmp.X_idx = Dvar_tmp.X_idx + 1
            Tvar_tmp.anchestor_x_idx = Tvar_tmp.anchestor_x_idx + 1
            Tvar_tmp.descendant_x_idx = Tvar_tmp.descendant_x_idx + 1
            
            Dvar_tmp['X'] = np.arange(1,len(Dvar_tmp)+1)
            
            sink = Dvar_tmp.X.max()+1
            
            # Inizialize graph
            G = nx.DiGraph()
            
            E = [] # Edges
            n = sink+1            
                
            for h in num_hyb:
                for idx,row in Dvar_tmp[(Dvar_tmp.hyb==h)].iterrows():
                    if h==1:
                        E.append((0,row.X, {'capacity': 1,'weight':0}))
                    # Add detection edges
                    E.append((row.X,n, {'capacity': 1,'weight':np.round(row.E_1*1000000).astype(int)}))
                    n = n+1
                    
                G.add_edges_from(E)
                E = []
                for idx, row in Tvar_tmp[(Tvar_tmp.anchestor_x_idx.isin(Dvar_tmp[Dvar_tmp.hyb==h].X_idx))].iterrows():
                    # Add transition edges
                    E.append((list(G.successors(row.anchestor_x_idx))[0],row.descendant_x_idx, {'capacity': 1,'weight':np.round(row.E_1*1000000).astype(int)}))
                G.add_edges_from(E)
                
            # For each D of last cycle connect to sink
            E = []
            for idx,row in Dvar_tmp[(Dvar_tmp.hyb==num_hyb.max())].iterrows():
                E.append((list(G.successors(row.X_idx))[0],sink, {'capacity': 1,'weight':0}))
            G.add_edges_from(E)
            
            # Prune graph removing leaf nodes
            remove_nodes = []
            for n in G.nodes:
                n_set = nx.algorithms.descendants(G, n)
                if not sink in n_set:
                    remove_nodes.append(n)
                    if n == 0: # source and sink are not connected
                        return
            
            remove_nodes.remove(sink)
            G.remove_nodes_from(remove_nodes)
            
            MaxFlowMinCost = nx.max_flow_min_cost(G, 0, sink)
            # Decode sequence
            E = []
            for n1 in MaxFlowMinCost:
                for n2 in MaxFlowMinCost[n1]:
                    if MaxFlowMinCost[n1][n2]==1:
                        E.append((n1,n2))
            G = nx.Graph()
            G.add_edges_from(E)
            G.remove_node(0)
            G.remove_node(sink)
            
            return {'G': G, 'Dvar': Dvar_tmp, 'Tvar': Tvar_tmp}
        

def runMaxFlowMinCost(data, d_th, k1, k2, tagList, n_threads, dth_max, prior):
    print("Generate Graph Model..."+" ("+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+") "+"\n")
    num_hyb=np.arange(1,int(np.amax(data.hyb))+1)
    data.sort_values('hyb', inplace=True)
    data=data.reset_index(drop=True)
    ## Graphical Model Data Structures
    print("[INFO] Add D vars"+" ("+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+") ")
    G = nx.Graph()
    G.add_nodes_from(data.index.values)    
    print("[INFO] Add T vars"+" ("+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+") ")
    for h1 in tqdm(num_hyb):
        KDTree_h1 = KDTree(data[data.hyb==h1][['x','y','z']])
        for h2 in num_hyb[h1:]:
            if h1 != h2:
                KDTree_h2 = KDTree(data[data.hyb==h2][['x','y','z']])
                query = KDTree_h1.query_ball_tree(KDTree_h2,d_th, p=2)    
                E = []
                offset1 = data.index[data.hyb==h1].min()
                offset2 = data.index[data.hyb==h2].min()
                for i,e1 in enumerate(query):
                    if e1:
                        for e2 in e1:
                            E.append((i+offset1,e2+offset2))
                G.add_edges_from(E)
        
    print("[INFO] Generate model"+" ("+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+") ")
    conn_comps = [list(c) for c in nx.connected_components(G)]
    for c in tqdm(range(len(conn_comps))):
        data.loc[conn_comps[c],'conn_comp']=c

    # Drop conn components with less than n_hybs elements
    gr = data.groupby('conn_comp')
    for i, group in gr:
         if len(group)<len(num_hyb):
             data=data.drop(group.index)    
    labels = np.unique(data.conn_comp)
    
    if labels.size>0:
        res = []
        for l in tqdm(np.nditer(labels),total=len(labels)):
            res.append(runComponent(data,int(l), d_th, k1, k2, tagList, num_hyb, dth_max, prior))
        #return maxFlowMinCost
        return [x for x in res if x is not None]