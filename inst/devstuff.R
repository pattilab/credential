#devstuff.R

if (F){
  setwd("X:/Nathaniel Mahieu/Q Exactive Plus/2NM107C Annotation Attempt")
  devtools::load_all("C:/Dropbox/GitHub/credential/")
  load("new_an.Rdata"); an;
  
  r.a = 0.9
  r.b = 1.1
  mzdiff = atm$c13 - atm$c12
  rt.lim = 10
  rcorr.lim = 0.7
  ppm.lim = 2.5
  mpc.f= 1.1
  charges = c(1,2,3,4,5,6) 
  
  pwms = credential(
    an,
    r.a = 0.9,
    r.b = 1.1,
    mzdiff = atm$c13 - atm$c12,
    rt.lim = 8,
    rcorr.lim = 0.7,
    ppm.lim = 2.5,
    mpc.f= 1.1,
    charges = c(1,2,3,4,5,6), 
  )
}
