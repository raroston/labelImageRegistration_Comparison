# labelImageRegistration_Comparison

## Prequisite:
ANTsR library and its prequisites. There are multiple ways to install ANTsR. Please follow their [**Installation** instructions on their repository](https://github.com/ANTsX/ANTsR?tab=readme-ov-file#installation). Additionally these R packages are necessary to run the scripts correctly.

`install.packages(c("stringr", "dplyr", "ggplot2", "ggsignif", "viridis"))`


## To generate results
All scripts are expected to be run inside the cloned repository using the syntax below. Note that when run, these scripts generate  ~40GB of output, so make sure there is enough space.

1. To run registration scripts: 
	* `Rscript ./Scripts/registration_intensityOnly.R`
	* `Rscript ./Scripts/registration_labelIntensity.R`
2. To run TBM analysis script:
	* `Rscript ./Scripts/analysis_TBM.R`
3. To run the PCA script:
	* `Rscript ./Scripts/analysis_PCA.R`
4. To run analysis of organ volumes:
	* Note: this script does not save plots to drive
	* `Rscript ./Scripts/analysis_organVolumes.R`
