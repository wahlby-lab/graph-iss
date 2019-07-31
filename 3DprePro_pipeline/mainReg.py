import numpy as np
import SimpleITK as sitk
from skimage import io
import pandas as pd
import sys
sys.path.insert(0, sys.argv[1])
from functions.parameters_map import SetRegParameters

#----------------------------------------------------------------------------------------------------------------------------------------------------- 
# READING FILES
#----------------------------------------------------------------------------------------------------------------------------------------------------- 
path = sys.argv[2]
path_csv = sys.argv[3] # csv file path
tile_csv=pd.read_csv(path_csv,sep='\t')
n_resolutions = str(sys.argv[4])
transform_type = str(sys.argv[5])
useNuclei = str(sys.argv[6])
if len(sys.argv)>=9:
    min_res = int(sys.argv[7])
    max_res = int(sys.argv[8])
    bspline_reg = str(sys.argv[9])
    n_resolutions = str(int(np.log2(min_res)-np.log2(max_res))+1)
    res_list = [] 
    tmp_res = min_res 
    while tmp_res>=max_res: 
        res_list.append(int(tmp_res)) 
        tmp_res = tmp_res/2
else:
    bspline_reg = str(sys.argv[7])

channels = 6
cycles = int(len(tile_csv)/channels)

# Find DO with z max and take it as reference
z_max = 0
for i in range(0,len(tile_csv),6):
    DO_tmp = np.squeeze(io.imread(tile_csv.File[i]))
    if DO_tmp.shape[0] > z_max:
        DO_ref = DO_tmp
        i_ref = i
        z_max = DO_tmp.shape[0]
del(DO_tmp)
io.imsave(path+'/reg_ch0'+str(int(i_ref/channels))+'.tif',DO_ref)

rows = DO_ref.shape[1]
cols = DO_ref.shape[2]
nbit=str(DO_ref.dtype)

if useNuclei =='Nuclei':
    DO_ref = sitk.GetImageFromArray(np.maximum(DO_ref,np.squeeze(io.imread(tile_csv.File[i_ref+1]))))
else:
    DO_ref = sitk.GetImageFromArray(DO_ref)


for i in range(0,len(tile_csv),6):
    cycle=int(i/channels)
    # Register curr_DO to DO_ref
    if useNuclei == 'Nuclei':
        curr_DO = sitk.GetImageFromArray(np.maximum(np.squeeze(io.imread(tile_csv.File[i])),np.squeeze(io.imread(tile_csv.File[i+1]))))
    else:
        curr_DO = sitk.GetImageFromArray(np.squeeze(io.imread(tile_csv.File[i])))

    if i!=i_ref:
        print("[INFO] DO"+str(cycle)+" registration...............................................................................................")
        elastixImageFilter2 = sitk.ElastixImageFilter()
        elastixImageFilter2.SetFixedImage(DO_ref)
        elastixImageFilter2.SetMovingImage(curr_DO)
        parameterMapVector = sitk.VectorOfParameterMap()
        for t in transform_type.split(","):
            parametersMap2=sitk.GetDefaultParameterMap(t)
            parametersMap2=SetRegParameters(parametersMap2)
            parametersMap2["NumberOfResolutions"]=[n_resolutions]
            if len(sys.argv)>=9:
                res_schedule = res_list * 3
                res_schedule.sort(reverse=True)
                res_schedule = [str(i) for i in res_schedule]
                parametersMap2["ImagePyramidSchedule"] = res_schedule
            parameterMapVector.append(parametersMap2)
    
        # Bspline reg
        if bspline_reg=="Bspline":
            parametersMap_deformable=sitk.GetDefaultParameterMap('bspline')
            parametersMap_deformable=SetRegParameters(parametersMap_deformable)
            parametersMap_deformable["FinalGridSpacingInVoxels"] = ['32', '32', '2']
            parametersMap_deformable["FinalGridSpacingInPhysicalUnits"]=""
            parametersMap_deformable["GridSpacingSchedule"]=""
            parametersMap_deformable["MovingImageDimension"]=['3']
            parametersMap_deformable["FixedImageDimension"]=['3']
            parametersMap_deformable["NumberOfResolutions"]=['1']
            if len(sys.argv)>=9:
                res_schedule = [res_list[-1]] * 3 
                res_schedule = [str(i) for i in res_schedule]
                parametersMap_deformable["ImagePyramidSchedule"] = res_schedule
            parametersMap_deformable["FinalBSplineInterpolationOrder"]=['0']
            parametersMap_deformable["MaximumNumberOfIterations"]= ['2000']
            parametersMap_deformable["UseCyclicTransform"]=['true']
            parameterMapVector.append(parametersMap_deformable)

        #elastixImageFilter2.SetParameterMap(parametersMap_deformable)
        transformixImageFilter2 = sitk.TransformixImageFilter()
        elastixImageFilter2.SetParameterMap(parameterMapVector)
        sitk.PrintParameterMap(parameterMapVector)
        elastixImageFilter2.Execute()
        transformParameterMap = elastixImageFilter2.GetTransformParameterMap()
        transformixImageFilter2.SetTransformParameterMap(transformParameterMap)
        del(elastixImageFilter2)    
        # Tranform and save registred curr_DO
        if useNuclei == 'Nuclei':
            curr_DO = sitk.GetImageFromArray(np.squeeze(io.imread(tile_csv.File[i])))

        print("[INFO] Transformix")
        transformixImageFilter2.SetMovingImage(curr_DO)
        transformixImageFilter2.Execute()
        reg_ch = sitk.GetArrayFromImage(transformixImageFilter2.GetResultImage())
        reg_ch[reg_ch<0] = 0
        reg_ch = reg_ch.astype(nbit)
        io.imsave(path+'/reg_ch0'+str(cycle)+'.tif',reg_ch)
        del(reg_ch)
        del(transformixImageFilter2)
    # registration (ch to respective DO)
    for ch in range(2,channels):
        # Do not register Nuclei
        if ch!=1:
            print("[INFO] DO"+str(cycle)+" ch"+str(ch)+" registration...............................................................................................")
            ch_toreg=sitk.GetImageFromArray(np.squeeze(io.imread(tile_csv.File[i+ch])))
            parametersMapIntraCycle = sitk.GetDefaultParameterMap('translation')
            parametersMapIntraCycle=SetRegParameters(parametersMapIntraCycle)
            parametersMapIntraCycle["NumberOfResolutions"]=['2']
            parametersMapIntraCycle["AutomaticParameterEstimation"]=["true"]
            reg_ch = sitk.GetArrayFromImage(sitk.Elastix(curr_DO,ch_toreg,parametersMapIntraCycle, True))  # Execute registration
            del(ch_toreg)
            if i == i_ref:
                reg_ch[reg_ch<0] = 0
                reg_ch = reg_ch.astype(nbit)
                io.imsave(path+'/reg_ch'+str(ch)+str(cycle)+'.tif',reg_ch)
                del(reg_ch)
            else:
                transformixImageFilter2 = sitk.TransformixImageFilter()
                transformixImageFilter2.SetTransformParameterMap(transformParameterMap)
                transformixImageFilter2.SetMovingImage(sitk.GetImageFromArray(reg_ch))
                transformixImageFilter2.Execute()
                reg_ch = sitk.GetArrayFromImage(transformixImageFilter2.GetResultImage())
                reg_ch[reg_ch<0] = 0
                reg_ch = reg_ch.astype(nbit)
                io.imsave(path+'/reg_ch'+str(ch)+str(cycle)+'.tif',reg_ch)
                del(reg_ch)
                del(transformixImageFilter2)
