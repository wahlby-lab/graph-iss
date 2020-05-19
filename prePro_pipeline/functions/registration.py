import numpy as np
import SimpleITK as sitk
from skimage import io
from functions.parameters_map import SetRegParameters
    
def RegistrationDO1(images,parametersMap,elastixImageFilter,transformixImageFilter,cycles,channels,path,save,nbit):
    """this registration is used when no general stain channel is available. Thus, MIP are used to achieve the alignment
        images=matrix with dimensions [cycles,channels,rows,cols] containing images to be aligned in stack
        parametersMap=registration type and parameters setting are set in the main script
        elastixImageFilter= elastix filter for registration
        transformixImageFilter= transformix filter to extract transformation matrix
        save= 0 -> don't save images
              1 -> save images"""
    
    reg_imgs=np.zeros(images.shape).astype(np.int32)  
    reg_imgs[0,:,:,:]=images[0,:,:,:]
    ref_DO=sitk.GetImageFromArray(images[0,0,:,:])
        
    elastixImageFilter.SetFixedImage(ref_DO)   
    
    for cycle in range(cycles):
        if all([save==1 , cycle == 0]):
            io.imsave(path+'/reg_ch00.tif',reg_imgs[0,0,:,:])
        reg_imgs[cycle,:,:,:]=images[0,:,:,:]
        for ch in range(2,channels):
            ch_toreg=sitk.GetImageFromArray(images[cycle,ch,:,:])
            transformixImageFilter.SetMovingImage(ch_toreg)
            parametersMap1 = parametersMap 
            parametersMap1['FinalBSplineInterpolationOrder'] = ['1']
            reg_ch = sitk.Elastix(ref_DO,ch_toreg,parametersMap1, True)
            reg_imgs[cycle,ch,:,:]=sitk.GetArrayFromImage(reg_ch)
            
            if save==1:
                io.imsave(path+'/reg_ch'+str(ch)+str(cycle)+'.tif',reg_imgs[cycle,ch,:,:])
            
    return reg_imgs        
            
def RegistrationDO(images,parametersMap,elastixImageFilter,transformixImageFilter,cycles,channels,path,save,nbit):
    """ this registration is used when all the general stain channel are available 
        images=matrix with dimensions [cycles,channels,rows,cols] containing images to be aligned in stack
        parametersMap=registration type and parameters setting are set in the main script
        elastixImageFilter= elastix filter for registration
        transformixImageFilter= transformix filter to extract transformation matrix
        save= 0 -> don't save images
              1 -> save images"""
    
    reg_imgs=np.zeros(images.shape).astype(nbit)         
    
    ref_DO=sitk.GetImageFromArray(images[0,0,:,:]) #first DO taken as reference and saved in the matrix 
    reg_imgs[0,0,:,:]=images[0,0,:,:] #saving DO first cycle (no MIP)
    
    #first registration (ch to respective DO)
    for cycle in range (cycles):
        if cycle==0:
            curr_DO=ref_DO
            if save==1:
                io.imsave(path+'/reg_ch00.tif',reg_imgs[0,0,:,:])
        else:
            curr_DO=sitk.GetImageFromArray(images[cycle,0,:,:]) #the first DO (cycle=0) is manatained unchanged (is the reference)
        for ch in range(2,channels):
            ch_toreg=sitk.GetImageFromArray(images[cycle,ch,:,:])
            parametersMapIntraCycle = parametersMap
            parametersMapIntraCycle["NumberOfResolutions"]=['2']
            reg_ch = sitk.Elastix(curr_DO,ch_toreg,parametersMapIntraCycle, True)
            reg_imgs[cycle,ch,:,:]=sitk.GetArrayFromImage(reg_ch)
            if all([save==1 , cycle == 0]): ######
                io.imsave(path+'/reg_ch'+str(ch)+'0.tif',reg_imgs[cycle,ch,:,:]) #to save channels registered images of the first cycle
    
    ref_DO=sitk.GetImageFromArray(images[0,0,:,:]) #first DO taken as reference and saved in the matrix
    elastixImageFilter.SetFixedImage(ref_DO)   
    for cycle in range(1,cycles):
        DO_toreg=sitk.GetImageFromArray(images[cycle,0,:,:]) #register DO
        elastixImageFilter2 = sitk.ElastixImageFilter()
        parametersMap2=sitk.GetDefaultParameterMap('rigid')
        parametersMap2=SetRegParameters(parametersMap2)
        parametersMap2["NumberOfResolutions"]=['2']
        elastixImageFilter2.SetParameterMap(parametersMap2)
        sitk.PrintParameterMap(parametersMap2)
        transformixImageFilter2 = sitk.TransformixImageFilter()
        elastixImageFilter2.SetFixedImage(ref_DO)
        elastixImageFilter2.SetMovingImage(DO_toreg)  #from the DO compared to the ref_DO we obtain transform
        elastixImageFilter2.Execute()
        
        transformParameterMap2 = elastixImageFilter2.GetTransformParameterMap()   
        transformixImageFilter2.SetTransformParameterMap(transformParameterMap2)       
            
        for ch in range(channels-1):
            if ch==0:
                ch_toreg=sitk.GetImageFromArray(images[cycle,0,:,:])
            else:
                ch+=1 #this is to avoid nuclei channel ch will be in [0,2,3,4,5]
                ch_toreg=sitk.GetImageFromArray(reg_imgs[cycle,ch,:,:])
                
            transformixImageFilter2.SetMovingImage(ch_toreg)
            transformixImageFilter2.Execute()
            reg_ch = transformixImageFilter2.GetResultImage()
            reg_imgs[cycle,ch,:,:]=sitk.GetArrayFromImage(reg_ch)
            if save==1:
                io.imsave(path+'/reg_ch'+str(ch)+str(cycle)+'.tif',reg_imgs[cycle,ch,:,:])
    return reg_imgs
