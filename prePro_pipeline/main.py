import numpy as np
import SimpleITK as sitk
from skimage import io
import pandas as pd
import  pickle
import tables
import sys
sys.path.insert(0, sys.argv[1])
from functions.parameters_map import SetRegParameters
from functions.registration import RegistrationDO, RegistrationDO1
from functions.candidates_extraction import ExtractCandidates
from functions.find_correspondences_create_inputs import findCorrespondencesMax_createInputs

#-----------------------------------------------------------------------------------------------------------------------------------------------------
#PARAMETERS SETTINGS
#-----------------------------------------------------------------------------------------------------------------------------------------------------
path = sys.argv[3] #output folder
path_csv=sys.argv[2] # csv file path
n_threads=int(sys.argv[6])
h_th=np.float(sys.argv[7])
sigma=np.float(sys.argv[8])
tile_csv=pd.read_csv(path_csv,sep='\t')
#-----------------------------------------------------------------------------------------------------------------------------------------------------
channels = 6
cycles = int(len(tile_csv)/channels)+1
final = []
for i in range(cycles):
    final.append(np.squeeze(io.imread(tile_csv.File[i]))) #DO
    final.append(np.zeros(np.squeeze(io.imread(tile_csv.File[i])).shape)) # Empty Nuclei
    final.append(np.squeeze(io.imread(tile_csv.File[i+cycles]))) #T
    final.append(np.squeeze(io.imread(tile_csv.File[i+cycles*2]))) #G
    final.append(np.squeeze(io.imread(tile_csv.File[i+cycles*3]))) #C
    final.append(np.squeeze(io.imread(tile_csv.File[i+cycles*4]))) #A

rows = final[0].shape[0]
cols = final[0].shape[1]
nbit=str(final[0].dtype)
radius=3

if np.amin([np.amax(final[x]) for x in range(0,channels*cycles,channels)])==0:
    sys.exit(0)

final=np.array(final).reshape(cycles,channels,rows,cols).astype(nbit) #contains all the channel and cycle of the same tile (cy,ch,rows,cols)
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#IMAGES REGISTRATION
#-----------------------------------------------------------------------------------------------------------------------------------------------------
try:
    elastixImageFilter = sitk.ElastixImageFilter()
    parametersMap=sitk.GetDefaultParameterMap('translation')
    parametersMap=SetRegParameters(parametersMap)
    elastixImageFilter.SetParameterMap(parametersMap)
    transformixImageFilter = sitk.TransformixImageFilter()
 
    if str(sys.argv[12])=="DO":
        final_reg=RegistrationDO(final,parametersMap,elastixImageFilter,transformixImageFilter,cycles,channels,path,1,nbit) # 0 - not save, 1 - save
    else:
        final_reg=RegistrationDO1(final,parametersMap,elastixImageFilter,transformixImageFilter,cycles,channels,path,1,nbit)
    final_reg=final_reg.astype(nbit)
except RuntimeError:
    print("ERROR: Registration")
    final_reg = final
    final_reg=final_reg.astype(nbit) 

#-----------------------------------------------------------------------------------------------------------------------------------------------------
#IMAGE NORMALIZATION
#-----------------------------------------------------------------------------------------------------------------------------------------------------
print("[INFO] Image Normalization....................................................................................................................")
mode_array = pickle.load(open(sys.argv[9], 'rb'))
#if mode array already exists
final_norm=np.zeros(final_reg.shape).astype(np.float64) #initialization final matrix of normalized images
for cycle in range(cycles):
    for ch in range(channels):
        if ch!=1:
            final_norm[cycle,ch,:,:]=((final_reg[cycle,ch,:,:].astype(np.int32)-mode_array[cycle,ch,0]).astype(np.float64))/(mode_array[cycle,ch,1] - mode_array[cycle,ch,0])
hdf5_path = sys.argv[10]
hdf5_file = tables.open_file(hdf5_path, mode='w')
data_storage = hdf5_file.create_array(hdf5_file.root, 'final_norm', final_norm)
hdf5_file.close()

final_norm[final_norm<0]=0
#CANDIDATES EXTRACTION
#-----------------------------------------------------------------------------------------------------------------------------------------------------
print("[INFO] Candidates Extraction....................................................................................................................")
candidates_max=np.zeros(final_norm.shape).astype('uint16')
#higher_contrast_points=np.zeros(final_norm.shape).astype('uint16') #mask with max points having an higher contrast (noise is discarded)

for cycle in range(cycles):
    for ch in range(channels-1):
        if ch>=1: #this is to avoid nuclei channel ch will be in [0,2,3,4,5]
            ch+=1
        if np.amax(final_norm[cycle,ch,:,:])==0:
            sys.exit(0)
        else:
            try:
                if ch==0:
                    candidates_max[cycle,ch,:,:]=ExtractCandidates(final_norm[cycle,ch,:,:], h_th*0.9, radius,nbit, sigma)
                    print(len(candidates_max[cycle,ch,:,:][candidates_max[cycle,ch,:,:]!=0]))
                else:
                    candidates_max[cycle,ch,:,:]=ExtractCandidates(final_norm[cycle,ch,:,:], h_th, radius,nbit, sigma)
                    print(str(cycle)+" "+str(ch)+" "+str(len(candidates_max[cycle,ch,:,:][candidates_max[cycle,ch,:,:]!=0])))
            except ValueError:
               print("ERROR: Extract candidates")
               sys.exit(0)

hdf5_path = sys.argv[11]
hdf5_file = tables.open_file(hdf5_path, mode='w')
data_storage = hdf5_file.create_array(hdf5_file.root, 'candidates_max', candidates_max)
hdf5_file.close()


#-----------------------------------------------------------------------------------------------------------------------------------------------------
#SIGNALS MERGING AND INPUT CREATION (IF NO CROSS-TALK COMPENSATION)
#-----------------------------------------------------------------------------------------------------------------------------------------------------
print("[INFO] Signal Merging..........................................................................................................................")
res=findCorrespondencesMax_createInputs(candidates_max,final_norm,cycles,channels,nbit,n_threads,radius)
inputs_df=res['inputs_df']
max_df=res['max_df']

inputs_df.to_pickle(sys.argv[4])
max_df.to_pickle(sys.argv[5])