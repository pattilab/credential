
#Credentialed Features

Credential metabolomic datasets based on isotopic signatures.

##Reference
http://pubs.acs.org/doi/abs/10.1021/ac503092d


##Installation

    install.packages("devtools")
    library(devtools)
    
    install_github("pattilab/credential")



##Example Usage

    library(credential)

    a = 1
    b = 2

    credentialed_features = credential(
      xs_a = loadXs(list.files(pattern=".xcmsSet$")[a]), # xcmsSet Object
      r_12t13_a = 10/9, # 12C/13C signal
      an_a = loadAn(list.files(pattern=".xsAnnotate$")[a]), # xsAnnotate Object
  
      xs_b = loadXs(list.files(pattern=".xcmsSet$")[b]), # xcmsSet Object
      r_12t13_b = 9/10, # 12C/13C signal
      an_b = loadAn(list.files(pattern=".xsAnnotate$")[b]), # xsAnnotate Object
  
      isotope_rt_delta_s = 5,
      ppm_for_isotopes = 5,
			ppm_for_isotopes_function=NULL,
  
  
      mixed_ratio_factor = 4,
      mixed_ratio_ratio_factor = 1.8,
      xs_a_file=list.files(pattern=".mzXML$")[a], # Path to raw data
      xs_b_file=list.files(pattern=".mzXML$")[b], # Path to raw data
      mpc_f = 1
    )
		
##Modifications Since Publication
 - The default, maximum mass-per-carbon limit was increased to accomidate heavier adducts such as Mg2+.  Previously these adducts could be improperly excluded if they increased the detected mass per carbon of a ner-limit compound.
 - The ppm error is now specifiable in a mass-dependent manner.  This allows the filtering to account for decreases in resolution at higher mass and the associated mass error introduced due to multiple peaks contributing to a single centroid mass.  For example the following function could be passed to ppm_for_isotopes_function:
		ppm_f = function(masses) {
			resolution = (2E6 * masses ^ -0.49) # Emperically determined mass dependence of resolution
			fwhm = masses / resolution
			ppm_plusminus = fwhm / 2 / masses * 1E6 + 0.3
			}