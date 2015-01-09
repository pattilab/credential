
#Credentialed Features

Credential metabolomic datasets based on isotopic signatures.

##Reference
http://pubs.acs.org/doi/abs/10.1021/ac503092d


##Installation
````r
install.packages("devtools")
library(devtools)

install_github("pattilab/credential")
````


##Example _credential()_ Usage
````r
library(credential)

# an accepts an xsAnnotate object from the CAMERA package.  This must contain grouping information between the two peaks as well. See below.

an = xsAnnotate(xs)

credentialed_features = credential(
	an = an,
	r_12t13_a = 0.97, # 12C/13C signal
	r_12t13_b = 1.23, # 12C/13C signal
	
	isotope_rt_delta_s = 5,
	ppm_for_isotopes = 5,

	mixed_ratio_factor = 4,
	mixed_ratio_ratio_factor = 1.8,
	mpc_f = 1.1,
	write_files=T
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
 - The default, maximum mass-per-carbon limit was increased to accomidate heavier adducts such as Mg2+.  Previously these adducts could be improperly excluded if they increased the detected mass per carbon of a ner-limit compound.
 - The ppm error is now specifiable in a mass-dependent manner.  This allows the filtering to account for decreases in resolution at higher mass and the associated mass error introduced due to multiple peaks contributing to a single centroid mass.  For example the following function could be passed to `ppm_for_isotopes_function`:
````r
ppm_f = function(masses) {
    resolution = (2E6 * masses ^ -0.49) # Emperically determined mass dependence of resolution
    fwhm = masses / resolution
    ppm_plusminus = fwhm / 2 / masses * 1E6 + 0.3
    }
````
 - Input has been simplified to better match the traditional xcms/CAMERA workflow.  Input is now a single xsAnnotate object containing two files.  This object must have been appropriately grouped and retention time corrected by the user as fitting the specific workflow.