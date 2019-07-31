#function to edit and set regitration map parameters 

def SetRegParameters(parametersMap):
    parametersMap["FixedImagePyramid"]= ["FixedRecursiveImagePyramid"]#
    parametersMap["ImageSampler"]=["RandomCoordinate"]
    parametersMap["UseRandomSampleRegion"]=["true"]
    parametersMap["MaximumNumberOfIterations"]= ['500']
    parametersMap["Metric"]= ["AdvancedNormalizedCorrelation"] 
    parametersMap["MovingImagePyramid "]= ["MovingRecursiveImagePyramid"]#
    parametersMap["NumberOfResolutions"]=['2'] 
    parametersMap["NumberOfSpatialSamples"]=['10000'] 
    parametersMap["NewSamplesEveryIteration"]=["true"] 
    parametersMap["BSplineInterpolationOrder"]=['1'] 
    parametersMap["FinalBSplineInterpolationOrder"]=['1']
    parametersMap["NumberOfHistogramBins"]=['128'] 
    parametersMap["Optimizer"]= ["AdaptiveStochasticGradientDescent"]
    parametersMap["ResampleInterpolator"]=["FinalBSplineInterpolator"]
    parametersMap["ResultImageFormat"]=["tif"]
    parametersMap["Interpolator"]=["BSplineInterpolator"]
    parametersMap["Registration"]=["MultiResolutionRegistration"]
    parametersMap["Resampler"]=["DefaultResampler"]
    parametersMap["HowToCombineTransforms "]=["Compose"]
    parametersMap["WriteResultImage"]= ["false"]
    parametersMap["FixedInternalImagePixelType"]=["long"]
    parametersMap["MovingInternalImagePixelType"]= ["long"]
    parametersMap["FixedImageDimension"]=['2']
    parametersMap["MovingImageDimension"]=['2']
    parametersMap["UseDirectionCosines"]=["true"]
    parametersMap["MinimumGradientMagnitude"]=['0.000001']
    parametersMap["MinimumStepLength"]=['0.01']
    parametersMap["AutomaticScalesEstimation"]=["true"]
    parametersMap["AutomaticParameterEstimation"]=["true"]
    parametersMap["AutomaticTransformInitialization"]=["true"]
    parametersMap["AutomaticTransformInitializationMethod"]=["GeometricalCenter"]
    parametersMap["HowToCombineTransforms"]=["Compose"]
    parametersMap["ErodeMask"]=["false"]
    parametersMap["ResultImagePixelType"]=["long"]
    return parametersMap

