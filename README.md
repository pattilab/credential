
#Credentialed Features

Credential metabolomic datasets based on isotopic signatures.


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
  
  
      mixed_ratio_factor = 4,
      mixed_ratio_ratio_factor = 1.8,
      xs_a_file=list.files(pattern=".mzXML$")[a], # Path to raw data
      xs_b_file=list.files(pattern=".mzXML$")[b], # Path to raw data
      mpc_f = 1
    )