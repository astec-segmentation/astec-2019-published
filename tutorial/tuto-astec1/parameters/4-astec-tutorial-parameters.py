PATH_EMBRYO = '.'

EN = '2019-Tutorial100'

begin = 0
end = 20
delta = 1
raw_delay = 0

## General parameters for segmentation propagation

astec_sigma1 = 0.6  		
astec_sigma2 = 0.15 		
astec_h_min_min = 4
astec_h_min_max = 18   		

## Glace Parameters (if astec_membrane_reconstruction_method is set to 1 or 2):
## membrane_renforcement

astec_sigma_membrane = 0.9
astec_sensitivity = 0.99  
astec_manual = False     	
astec_manual_sigma = 15   
astec_hard_thresholding = False 
astec_hard_threshold = 1.0      

## Tensor voting framework

astec_sigma_TV = 3.6    
astec_sigma_LF = 0.9    
astec_sample = 0.2      

## Default parameters (for classical use, default values should not be changed)

astec_RadiusOpening = 20 		
astec_Thau = 25 				
astec_MinVolume = 1000 		
astec_VolumeRatioBigger = 0.5 
astec_VolumeRatioSmaller = 0.1
astec_MorphosnakeIterations = 10 
astec_NIterations = 200 		
astec_DeltaVoxels = 10**3  	
astec_nb_proc = 10
