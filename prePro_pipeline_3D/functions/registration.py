import numpy as np
import SimpleITK as sitk
from skimage import io
from functions.parameters_map import SetRegParameters
            
def RegistrationDO(final,parametersMap,elastixImageFilter,transformixImageFilter,cycles,channels,path,save,nbit):
    """ this registration is used when all the general stain channel are available 
        images=matrix with dimensions [cycles,channels,rows,cols] containing images to be aligned in stack
        parametersMap=registration type and parameters setting are set in the main script
        elastixImageFilter= elastix filter for registration
        transformixImageFilter= transformix filter to extract transformation matrix
        save= 0 -> don't save images
              1 -> save images"""
    z_max= max(map(len, final))    
    reg_imgs=np.zeros((cycles,channels,z_max,final[0].shape[1],final[0].shape[2])).astype(nbit)
    
    ref_DO=sitk.GetImageFromArray(final[0]) #first DO taken as reference and saved in the matrix 
    reg_imgs[0,0,0:final[0].shape[0],:,:]=final[0]
    
    for i in range(0,len(final),channels):
        cycle=int(i/channels)
        #first registration (ch to respective DO)
        if cycle==0:
            curr_DO=ref_DO
            if save==1:
                io.imsave(path+'/reg_ch00.tif',reg_imgs[0,0,:,:,:])
        else:
            curr_DO=sitk.GetImageFromArray(final[i])
        for ch in range(2,channels):
            print("[INFO] DO"+str(cycle)+" ch"+str(ch)+" registration...............................................................................................")
            ch_toreg=sitk.GetImageFromArray(final[i+ch])
            parametersMapIntraCycle = parametersMap
            parametersMapIntraCycle["NumberOfResolutions"]=['2']
            parametersMapIntraCycle["AutomaticParameterEstimation"]=["true"]
            reg_ch = sitk.Elastix(curr_DO,ch_toreg,parametersMapIntraCycle, True)
            reg_imgs[cycle,ch,0:sitk.GetArrayFromImage(reg_ch).shape[0],:,:]=sitk.GetArrayFromImage(reg_ch)
            if all([save==1 , cycle == 0]):
                io.imsave(path+'/reg_ch'+str(ch)+'0.tif',reg_imgs[cycle,ch,:,:,:]) #to save channels registered images of the first cycle
   
    for i in range(channels,len(final),channels):
        cycle=int(i/channels)
        print("[INFO] DO"+str(cycle)+" registration.................................................................................................................")
        ref_DO=sitk.GetImageFromArray(final[0])
        DO_toreg=sitk.GetImageFromArray(final[i]) #register DO+Nuclei 
        elastixImageFilter2 = sitk.ElastixImageFilter()
        elastixImageFilter2.SetFixedImage(ref_DO)
        elastixImageFilter2.SetMovingImage(DO_toreg)
        parametersMap2=sitk.GetDefaultParameterMap("rigid")
        parametersMap2=SetRegParameters(parametersMap2)
        parametersMap2["NumberOfResolutions"]=["2"]
        parametersMap_deformable=sitk.GetDefaultParameterMap('bspline')      
        parametersMap_deformable=SetRegParameters(parametersMap_deformable)   
        parametersMap_deformable["FinalGridSpacingInVoxels"] = ['32', '32', '2']
        parametersMap_deformable["FinalGridSpacingInPhysicalUnits"]=""
        parametersMap_deformable["GridSpacingSchedule"]=""
        parametersMap_deformable["MovingImageDimension"]=['3']
        parametersMap_deformable["FixedImageDimension"]=['3']
        parametersMap_deformable["NumberOfResolutions"]=['1']
        parametersMap_deformable["FinalBSplineInterpolationOrder"]=['0']
        parametersMap_deformable["MaximumNumberOfIterations"]= ['2000']
        parametersMap_deformable["UseCyclicTransform"]=['true']               
        elastixImageFilter2.SetParameterMap(parametersMap_deformable)
        sitk.PrintParameterMap(parametersMap_deformable)             
        transformixImageFilter2 = sitk.TransformixImageFilter()
        parameterMapVector = sitk.VectorOfParameterMap()
        parameterMapVector.append(parametersMap2)
        parameterMapVector.append(parametersMap_deformable)
        elastixImageFilter2.SetParameterMap(parameterMapVector)
        elastixImageFilter2.Execute()
        transformParameterMap = elastixImageFilter2.GetTransformParameterMap()
        transformixImageFilter2.SetTransformParameterMap(transformParameterMap)

        for ch in range(channels-1):
            if ch==0:
                ch_toreg=sitk.GetImageFromArray(final[i])
            else:
                ch+=1 #this is to avoid nuclei channel ch will be in [0,2,3,4,5]
                ch_tmp = reg_imgs[cycle,ch,:,:,:]
                ch_toreg=sitk.GetImageFromArray(ch_tmp[~(ch_tmp == 0).all(axis=(1,2))])
                
            transformixImageFilter2.SetMovingImage(ch_toreg)
            transformixImageFilter2.Execute()
            reg_ch = transformixImageFilter2.GetResultImage()
            reg_imgs[cycle,ch,0:sitk.GetArrayFromImage(reg_ch).shape[0],:,:]=sitk.GetArrayFromImage(reg_ch)

    for cycle in range(cycles):
        for z in range(z_max):
            if np.amax(reg_imgs[cycle,0,z,:,:])==0:
                reg_imgs[cycle,:,z,:,:]=np.zeros((reg_imgs[cycle,:,z,:,:].shape))

        if save==1:
            for ch in range(channels):
                io.imsave(path+'/reg_ch'+str(ch)+str(cycle)+'.tif',reg_imgs[cycle,ch,:,:,:])

    return reg_imgs
