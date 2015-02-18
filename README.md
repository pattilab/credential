#Credential
Credential metabolomic datasets based on isotopic signatures.

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
  mzdiff = atm$c13 - atm$c12, # The difference in mass between the enriched and natural abundance species
  rt.lim = 10, # A rough filter of RT. Make as large as possible
  rcorr.lim = 0.7, # The minimum pearson correlation between isotope EICs
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
  - This is version 2.0 of the credentialing software.  The principle of the process is identical but major routines have been rewritten for speed and breadth.

##Modifying The Source
 1. Clone this repo
 2. Modify the source
 3. Run:
 ````r
 devtools::load_all("path/to/credential/")
 ````