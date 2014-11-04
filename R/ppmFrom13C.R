ppm_from_13C = function(m1, m2, charge=1) {
  
  if (m2 < m1) {
    tmp = m2
    m2 = m1
    m1=tmp
  }
  
  isotope_predicted_carbons = round(abs(m1 - m2) * charge / (aC13-aC12))
  isotope_predicted_mass = isotope_predicted_carbons*(aC13-aC12)/charge + m1
  
  ppm_error = (m2 - isotope_predicted_mass) / isotope_predicted_mass * 1E6
  return(ppm_error)
}