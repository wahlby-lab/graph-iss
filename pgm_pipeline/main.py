import numpy as np
import pandas as pd
import itertools
import sys
sys.path.insert(0, sys.argv[1])
import pickle
import tables
from keras.models import model_from_json
from tqdm import tqdm

from basecalling.pgm_basecalling import basecallingMaxFlowMinCost, signalDecoding

try:
    dataset = pd.read_pickle(sys.argv[2])
    max_df = pd.read_pickle(sys.argv[11])
    read_hdf5_file = tables.open_file(sys.argv[12], mode='r')
    final_norm = read_hdf5_file.root.final_norm[:]
    read_hdf5_file.close()
    if len(dataset.columns)==6: # 2D dataset
        dataset['z']=np.zeros(len(dataset))
        dataset = dataset[['cycle', 'ch', 'x', 'y', 'z', 'Intensities_window_5x5', 'prob_DNN']]
        
        max_df['z_T']=np.zeros(len(max_df))
        max_df['z_G']=np.zeros(len(max_df))
        max_df['z_C']=np.zeros(len(max_df))
        max_df['z_A']=np.zeros(len(max_df))
        max_df = max_df[['I_T', 'I_G', 'I_C', 'I_A', 'x_T', 'y_T', 'z_T', 'x_G', 'y_G', 'z_G', 'x_C', 'y_C', 'z_C', 'x_A', 'y_A', 'z_A', 'cycle']]
        
        final_norm=final_norm[:,:,np.newaxis,:,:]
    dataset.drop_duplicates(['cycle','ch','x','y','z'], inplace=True)
except EOFError:
    f = open(sys.argv[4],'w')
    f.write('letters,global_X_pos,global_Y_pos,tile_ID,max_dist,seq_quality_min')
    f.close()
    sys.exit(0)

if sys.argv[9]=="prior":
    tags = pd.read_csv(sys.argv[10], sep = ",", usecols = [0], header = None, names = ["barcode_bases"])
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('A','5');
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('C','4');
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('G','3');
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('T','2');
    tags = tags['barcode_bases'].apply(lambda x: pd.Series(list(x)))
else:
    # Generate False Barcodes
    num_hybs = np.amax(dataset.cycle)+1
    tags = pd.DataFrame(list(itertools.product("ACGT", repeat=num_hybs)))
    tags = pd.DataFrame(tags.apply(lambda x: ''.join(x), axis=1), columns=['barcode_bases'])
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('A','5');
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('C','4');
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('G','3');
    tags['barcode_bases'] = tags['barcode_bases'].str.replace('T','2');
    tags = tags['barcode_bases'].apply(lambda x: pd.Series(list(x)))

## 4. Run basecalling (NN + pgm)
net_path=sys.argv[13]
model_name='/model_FT3'
json_file = open(net_path+model_name+'.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
DNN_model = model_from_json(loaded_model_json)
DNN_model.load_weights(net_path+model_name+'_weights')

d_th=float(sys.argv[6])
cycles=dataset.cycle.max()+1
try:
    res = signalDecoding(dataset, d_th, k1=0.33, k2=0.1, k3=0.4, k4=0.8, tagList=tags, n_threads = int(sys.argv[5]), dth_max=float(sys.argv[7]), prior=sys.argv[9])
except ValueError:
    sys.exc_info()
    sys.exit(0)

final_norm[final_norm<0]=0

print("Generating output file...")
barcodes_list=[]
signals_df_list=[]
if res is not None:
    for i in tqdm(range(len(res))):
        tmp_res=res[i]
        baseCalling_res = basecallingMaxFlowMinCost(tmp_res['G'], tmp_res['Dvar'], tmp_res['Tvar'], max_df, final_norm, DNN_model, num_hybs)
        barcodes = baseCalling_res['barcodes']
        signals_df = baseCalling_res['signals_df'] 
        T_Costs = baseCalling_res['T_Costs']
        D_Costs = baseCalling_res['D_Costs']
        for j in range(len(barcodes)):
            barcodes_list.append(np.hstack((barcodes.values[j],D_Costs[j],T_Costs[j])))
            signals_df_list.append(signals_df[j])
        
if len(barcodes_list):        
    barcodes_df = pd.DataFrame(barcodes_list)
    barcodes_df['barcode'] = barcodes_df.iloc[:,3:-2].apply(lambda x: ''.join(x), axis=1)
    barcodes_df = barcodes_df.drop(barcodes_df.columns[3:-3], axis=1)
    barcodes_df = barcodes_df.rename(index=str, columns={barcodes_df.columns[0]:"x",barcodes_df.columns[1]:'y',barcodes_df.columns[2]:'z',barcodes_df.columns[3]:'D_Cost',barcodes_df.columns[4]:'T_Cost'})
    barcodes_df.to_pickle(sys.argv[3],protocol=2)
    pickle.dump(signals_df_list, open(sys.argv[8], 'wb'))

    barcodes_csv=pd.DataFrame(barcodes_df[["barcode"]],copy=True);barcodes_csv['global_X_pos']=barcodes_df[['x']].values;barcodes_csv["global_Y_pos"]=barcodes_df[["y"]];barcodes_csv["tile_ID"]=np.ones((len(barcodes_df),1));barcodes_csv["max_dist"]=barcodes_df.T_Cost;barcodes_csv["seq_quality_min"]=barcodes_df.D_Cost;
    barcodes_csv=barcodes_csv.rename(index=str, columns={"barcode": "letters"})
    barcodes_csv.to_csv(sys.argv[4],index=False)