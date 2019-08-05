""" Course registration.

    Parameters
    ----------
    sys.argv[1] : Path to pre-processing library
    sys.argv[2] : Output folder where to store aligned images
    sys.argv[3] : Input csv file contaning an array of input images
    sys.argv[4] : Number of resolution levels used for registration
    sys.argv[5] : Type of registration (valid arguments: "translati-
        on", "rigid", or "affine")
    sys.argv[6] : String flag that enables sequencing rounds registra-
        tion using maximum projected images of general stain and
        nuclei channel if enabled. Genaral stain images only are used
        if disabled. (valid arguments: "Nuclei": feature enabled,
        other string: feature disabled)
    sys.argv[7] : Type of registration procedure. (valid arguments:
        "DO1" : if only the general stain of the first sequencing round
        is available
        "DO" : if a general stain image is available for each sequencing
        round
    sys.argv[8] : highest resolution level for multiresolution image reg-
        istration
    sys.argv[9] : lower resolution level for multiresolution image regi-
        stration
    sys.argv[10] : flag to enable BSpline registration after rigid
"""

import SimpleITK as sitk
from skimage import io
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, sys.argv[1])
from functions.parameters_map import SetRegParameters

#-----------------------------------------------------------------------------------------------------------------------------------------------------
#PARAMETERS SETTINGS
#-----------------------------------------------------------------------------------------------------------------------------------------------------
path = sys.argv[2] #output folder
path_csv=sys.argv[3] # csv file path
tile_csv=pd.read_csv(path_csv,sep='\t')
n_resolutions = str(sys.argv[4])
transform_type = str(sys.argv[5])
useNuclei = str(sys.argv[6])
reg_type = str(sys.argv[7])
if len(sys.argv)>9:
    min_res = int(sys.argv[8])
    max_res = int(sys.argv[9])
    bspline_reg = str(sys.argv[10])
    n_resolutions = str(int(np.log2(min_res)-np.log2(max_res))+1)
    res_list = []
    tmp_res = min_res
    while tmp_res>=max_res:
        res_list.append(int(tmp_res))
        tmp_res = tmp_res/2
else:
    bspline_reg= str(sys.argv[8])

channels = 6
cycles = int(len(tile_csv)/channels)

# Take DO0 as reference
DO_ref = np.squeeze(io.imread(tile_csv.File[0]))
io.imsave(path+'/reg_ch00.tif',DO_ref)
rows = DO_ref.shape[0]
cols = DO_ref.shape[1]
nbit=str(DO_ref.dtype)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# REGISTRATION
#----------------------------------------------------------------------------------------------------------------------------------------------------
if useNuclei =='Nuclei':
    DO_ref = sitk.GetImageFromArray(np.maximum(DO_ref,np.squeeze(io.imread(tile_csv.File[1]))))
else:
    DO_ref = sitk.GetImageFromArray(DO_ref)

if reg_type == 'DO1':
    elastixImageFilter = sitk.ElastixImageFilter()
    transformixImageFilter = sitk.TransformixImageFilter()
    parametersMap=sitk.GetDefaultParameterMap(transform_type)
    parametersMap=SetRegParameters(parametersMap)
    parametersMap['FinalBSplineInterpolationOrder'] = ['1']
    parametersMap["NumberOfResolutions"]=[n_resolutions]
    sitk.PrintParameterMap(parametersMap)
    elastixImageFilter.SetParameterMap(parametersMap)
    elastixImageFilter.SetFixedImage(DO_ref)
    #Register everything to DO0
    for i in range(0,len(tile_csv),6):
        cycle=int(i/channels)
        for ch in range(channels):
            if ch!=1:
                # Register curr_DO to DO_ref
               if useNuclei == 'Nuclei':
                   mov_ch = sitk.GetImageFromArray(np.maximum(np.squeeze(io.imread(tile_csv.File[i])),np.squeeze(io.imread(tile_csv.File[i+1]))))
               else:
                   mov_ch = sitk.GetImageFromArray(np.squeeze(io.imread(tile_csv.File[i])))
               elastixImageFilter.SetMovingImage(mov_ch)
               elastixImageFilter.Execute()
               transformParameterMap = elastixImageFilter.GetTransformParameterMap()
               transformixImageFilter.SetTransformParameterMap(transformParameterMap)
               if useNuclei == 'Nuclei':
                   mov_ch = sitk.GetImageFromArray(np.squeeze(io.imread(tile_csv.File[i])))
               transformixImageFilter.SetMovingImage(mov_ch)
               transformixImageFilter.Execute()
               reg_ch = sitk.GetArrayFromImage(transformixImageFilter.GetResultImage())
               reg_ch[reg_ch<0] = 0
               reg_ch = reg_ch.astype(nbit)
               io.imsave(path+'/reg_ch'+str(ch)+str(cycle)+'.tif',reg_ch)

else:
    for i in range(0,len(tile_csv),6):
        cycle=int(i/channels)
        # Register curr_DO to DO_ref
        if useNuclei == 'Nuclei':
            curr_DO = sitk.GetImageFromArray(np.maximum(np.squeeze(io.imread(tile_csv.File[i])),np.squeeze(io.imread(tile_csv.File[i+1]))))
        else:
            curr_DO = sitk.GetImageFromArray(np.squeeze(io.imread(tile_csv.File[i])))
        
        if i !=0:
            print("[INFO] DO"+str(cycle)+" registration...............................................................................................")
            elastixImageFilter2 = sitk.ElastixImageFilter()
            transformixImageFilter2 = sitk.TransformixImageFilter()
            elastixImageFilter2.SetFixedImage(DO_ref)
            elastixImageFilter2.SetMovingImage(curr_DO)
            parameterMapVector = sitk.VectorOfParameterMap()
            for t in transform_type.split(","):
                parametersMap2=sitk.GetDefaultParameterMap(t)
                parametersMap2=SetRegParameters(parametersMap2)
                parametersMap2['FinalBSplineInterpolationOrder'] = ['1']
                parametersMap2["NumberOfResolutions"]=[n_resolutions]
                if len(sys.argv)>=9:
                    res_schedule = res_list * 2
                    res_schedule.sort(reverse=True)
                    res_schedule = [str(i) for i in res_schedule]
                    parametersMap2["ImagePyramidSchedule"] = res_schedule
                parameterMapVector.append(parametersMap2)

            # BSpline reg
            if bspline_reg=="BSpline":
                parametersMap_deformable=sitk.GetDefaultParameterMap('bspline')
                parametersMap_deformable=SetRegParameters(parametersMap_deformable)
                parametersMap_deformable["FinalGridSpacingInVoxels"] = ['32', '32']
                parametersMap_deformable["FinalGridSpacingInPhysicalUnits"]=""
                parametersMap_deformable["GridSpacingSchedule"]=""
                parametersMap_deformable["MovingImageDimension"]=['2']
                parametersMap_deformable["FixedImageDimension"]=['2']
                parametersMap_deformable["NumberOfResolutions"]=['1']
                parametersMap_deformable["FinalBSplineInterpolationOrder"]=['1']
                parametersMap_deformable["MaximumNumberOfIterations"]= ['2000']
                if len(sys.argv)>=9:
                    res_schedule = [res_list[-1]] * 2
                    res_schedule = [str(i) for i in res_schedule]
                    parametersMap_deformable["ImagePyramidSchedule"] = res_schedule
                parameterMapVector.append(parametersMap_deformable)

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
                parametersMapIntraCycle=sitk.GetDefaultParameterMap('translation')
                parametersMapIntraCycle=SetRegParameters(parametersMapIntraCycle)
                parametersMapIntraCycle['FinalBSplineInterpolationOrder'] = ['1']
                parametersMapIntraCycle["NumberOfResolutions"]=['2']
                parametersMapIntraCycle["AutomaticParameterEstimation"]=["true"]
                reg_ch = sitk.GetArrayFromImage(sitk.Elastix(curr_DO,ch_toreg,parametersMapIntraCycle, True))  # Execute registration
                del(ch_toreg)
            if cycle == 0:
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
