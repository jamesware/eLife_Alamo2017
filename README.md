# Effect of Myosin Variants on Interacting-Heads Motif Explain Distinct Hypertrophic and Dilated Cardiomyopathy Phenotypes
### _Lorenzo Alamo\*, James S. Ware\*, Antonio Pinto, Richard E. Gillilan, Jonathan G. Seidman, Christine E. Seidman\* and Raúl Padrón\*_

_\* equal contributions_

This repository contains code to reproduce analyses presented in the paper above, and to generate the figures.  This code is released under a GPL3 license.

The whole directory is formatted as an RStudio project, with the project definition file `eLife_Alamo2017.Rproj`.  

The directory `Chimera` contains chimera scripts to reproduce structural images. The myosin interacting heads structure described in this paper is deposited in the Protein Data Bank with accession **5TBY**.  

The folder `R` contains code to reproduce statistical analyses.  To recapitulate these, open the project, and run `R/analysis/analysis.Rmd`.  This rmarkdown document will produce an html file containing all of the results, and will write the supplementary files as word documents in the `analysis` subfolder.  

The analysis script will also download the published data used in this analysis (using wget, called by a `downloadData()` function defined here).  These will be downloaded to `/data-raw`. This folder also contains a slice of the ExAC data required for these analyses (data saved to the repo for convenience).
