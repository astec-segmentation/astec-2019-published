# ACE  for Automated Cell Extractor : THIS FILE IS A PATTERN FOR THE USAGE OF GACE METHOD, THIS IS NOT SUPPOSED TO BE USED AS A STEP FOR THE ASTEC METHOD !! -> THE GACE METHOD USED FOR MARS ALGORITHM IS ALREADY WRITTEN IN MARS FILES...
from definitions import *
if not os.path.isdir(reconstruct_Path):
    os.mkdir(reconstruct_Path)  
os.system("cp -f "+astec_Path+"definitions.py "+reconstruct_Path )
os.system("cp -f "+astec_Path+"0-ace.py "+reconstruct_Path ) 

#Import Functions
from ImageHandling import imsave #, imread, SpatialImage
from lineage import timeNamed
from ACE import GACE

### Parameters:

# membrane_renforcement
sigma_membrane=0.9 # parametre de rehaussement de membranes (en unites reelles, a priori 0.9 um est adapte a des images fusionnees de Patrick/Ralph/Aquila)

# anisotropicHist /!\ etape critique
sensitivity=0.99 # parametre de binarisation des membranes, etape critique /!\ en cas d'echec, privilegier une parametrisation "manuelle" de la fonction anisotropicHist via l'activation de l'option 'manual'

manual=False     # par defaut, garder ce parametre a False. Si echec (ie seuils aberrants, se traduisant par une image binarisee des membranes de qualite tres mediocre), 
                 # le mettre a True et relancer les calculs sur l'image a tester. Si nouvel echec, jouer sur la valeur de manual_sigma... resultat non garanti
manual_sigma=15  # parametre d'ajustement initial des histogrammes axiaux permettant le calcul de seuils de binarisation de l'image des membranes (parametre pris en compte seulement si manual = True). 
                 # Faire varier manual_sigma entre 5 et 25 en cas d'echec au premier essai. Resultat non garanti...

hard_thresholding=False  # Si echec de la precedente methode de seuillage des membranes par calcul automatique de seuils directionnels, possibilite de choisir un seuillage dur global sur l'image en mettant cette option a True
hard_threshold=1.0       # Si hard_thresholding est a True, seuillage des membranes rehaussees via ce seuil (valeur 1 : adaptee pour le time-point t001 d'Aquila par exemple)

# TVmembrane
sigma_TV=3.6     # parametre definissant l'echelle de vote (ie de propagation) des membranes par vote de tenseurs (en reelles), a priori a choisir entre 3 um (cellules etroites) et 4.5 um (gros gaps a combler)
sigma_LF=0.9     # parametre de lissage de l'image reconstruite (en reelles), normalement la valeur par defaut = 0.9 um est satisfaisante
sample=0.2       # parametre permettant d'optimiser la vitesse de calcul du vote de tenseurs (eviter de toucher)



### PROCESS GLOBAL IMAGE RECONSTRUCTION VIA GACE METHOD (Global Automatic Cell Extractor)

for t in range(begin, end):
    time_reconstruction=t #Time point of Reconstruction
    print 'Starting the reconstruction at ' + str(time_reconstruction)
    fused_file=timeNamed(fused_files,time_reconstruction) #Actual image file to reconstruct
    reconstruct_file=timeNamed(reconstruct_files,time_reconstruction) #Output Reconstruction file
    
    # Global Automated Cell Extractor

    reconstruct_image=GACE(fused_file, sigma_membrane=sigma_membrane, sensitivity=sensitivity, manual=manual, manual_sigma=manual_sigma, hard_thresholding=hard_thresholding, hard_threshold=hard_threshold,
            sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample, keep_all=False, keep_membrane=False, verbose=True)

    #SAVE OUTPUT
    print 'Write the image reconstruction in ' + reconstruct_file
    imsave(reconstruct_file, reconstruct_image)

print 'GACE RECONSTRUCTION DONE'
