library(credential)

a = 1
b = 2

credentialed_features = credential(
  xs_a = loadXs(list.files(pattern=".xcmsSet$")[a]), 
  r_12t13_a = 10/9,
  an_a = loadAn(list.files(pattern=".xsAnnotate$")[a]), 
  
  xs_b = loadXs(list.files(pattern=".xcmsSet$")[b]), 
  r_12t13_b = 9/10,
  an_b = loadAn(list.files(pattern=".xsAnnotate$")[b]), 
  
  isotope_rt_delta_s = 5,
  ppm_for_isotopes = 5,
  
  
  mixed_ratio_factor = 4,
  mixed_ratio_ratio_factor = 1.8,
  xs_a_file=list.files(pattern=".mzXML$")[a], 
  xs_b_file=list.files(pattern=".mzXML$")[b],
  mpc_f = 1
)
