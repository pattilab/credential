
#Credentialed Features
#Credential metabolomic datasets based on isotopic signatures.

##Reference
http://pubs.acs.org/doi/abs/10.1021/ac503092d


##Install
````r
# install.packages("devtools")
devtools::install_github("pattilab/credential")
````


##Example _credential()_ Usage
````r
library(credential)

# an accepts an xsAnnotate object from the CAMERA package.  This an should have two experimental runs of credentialed standard extract.

credentialed_features = credential(
  an = an,
  r.a = 0.9, # Ratio of sample 1
  r.b = 1.1, # Ratio of sample 2
  mzdiff = aC13-aC12, # The difference in mass between the enriched and natural abundance species
  rt.lim = 10, # A rough filter of RT
  rcorr.lim = 0.7, # The minimum pearson correlation between isotopes
  ppm.lim = 2.5, # The ppm limit of isotopes
  mpc.f= 1.1, # A factor by which to expand the mass per carbon limits
  charges = c(1,2,3,4,5,6) # The charge states to consider
)
````

##Example Data Preparation
````r
library(xcms)
library(CAMERA)

xs = group(retcor(xcmsSet("file_a.mzxml", "file_b.mzxml")))
an = xsAnnotate(xs)

# settings for xcmsSet(), group(), and retcor() must be optimized by the user
# an is then passed to the credential() call as above
````

##Modifications Since Publication
  - Much has changed.
 
 ##Modifying The Source
 1. Clone this repo
 2. Modify the source
 ````r
 devtools::load_all("path/to/credential/")
 ````