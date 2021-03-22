## Deciphering how early life adiposity influences breast cancer risk using Mendelian randomization

_This repository contains materials for the mini-project I've done as a part of my PhD programme at the University of Bristol_

_Publication is curently in preparation_

This README contains:
- e-poster summarising the project research (2 Nov 2020)
- project background and aims
- outline of the project work/code stored in this repo


### E-poster created for NCRI Virtual conference 2020

<img src="e-poster_NCRI.png" width="1000"> 



### Project Background 

Multiple observational studies and a recent MR study (Richardson et al., 2020) have shown that early life adiposity may be protective against later-life breast cancer. The biological mechanism underlying this effect is unclear but is likely to independent of adiposity in adulthood. In this project, the aim is to identify potential mediators of the protective effect of high childhood BMI on breast cancer risk. 

As a part of this project I have performed the following (depending on data availability):

  -	Extensive literature study into potential mediators that may be influenced by childhood BMI and affect breast cancer risk ([lit review summary here][1])
  -	MR of childhood BMI on the identified potential mediators (step 1 of 2-step MR)
  -	MR of the identified potential mediators on breast cancer (step 2 of 2-step MR)
  -	Multivariable MR analysis on mediators that show causal effects in the 2-step MR
  -	Mediation analysis using the product and difference methods to calculate mediated direct and indirect effects of the individual mediators on breast cancer risk

The five main categories of potential mediators considered ([mediators summary table][2]):

  1.	Hormones (IGF1, oestradiol, SHBG (pre/post-menopause), testosterone (pre/post menopause; free/total/bioavailable))
  2.	Physical traits (mammographic density, breast size)
  3.	Reproductive traits (age at menarche, age at first birth, number of births
  4.  Glycaemic traits (fasting insulin, fasting glucose, HOMA-IR, HOMA-B, Hba1c)



I have created a script that automatically runs all MR analyses for a specified mediator and produces a report with all results and sensitivity tests/plots. [This folder][3] contains PDFs for all hormone mediators, with all sensitivity tests, and the summary of all MR analyses is presented on the last page. 



### This repository

Main analysis scripts and metadata (see details below):

```
├── set_paths.R
├── 01_process_gwas_summary.Rmd
├── 02_mr_BMI_to_mediators.Rmd
├── 03_mr_mediators-to-BC.Rmd
├── 04_mvmr_run_analysis.Rmd
├── 05_mvmr_create_plots.Rmd
├── 06_mediation_analysis.Rmd
├── 07_generate_MR_reports.Rmd
├── MRreport.Rmd
├── functions.R
├── functions_mvmr.R
├── metadata
│   └── data_lookup.csv
│   └── pheno_correlations.csv

```

Supplementary scripts:

```
├── code_supplementary
│   ├── extract_from_GWAScatalog.Rmd
│   └── review_mediators_in_MRBase.Rmd
│   └── combine_data_from_reports.Rmd
│   └── create_pheno_covariance_matrix.Rmd

```

Old scripts that were retired after tidy main analysis scripts were created:

```
├── code_legacy
│   ├── mvmr.R
│   ├── mvmr_combined_plots.Rmd
│   ├── performMR_BMI-to-BC.Rmd
│   ├── performMR_BMI-to-Mediators.Rmd
│   └── performMR_Mediators-to-BC.Rmd
```

Some exploratory scripts that may or may not be used for anything:

```
├── code_exploratory
│   ├── examine_phenotype_files.R
│   ├── experimental_MVMR.Rmd
│   ├── extract_instrumets_test.R
│   └── making_DAGs.R
```



### How to reproduce this analysis, i.e. Main Workflow


1.  Script `set_paths.R` is imported by all other scripts, and is used to set the environment of where the project is run.

2. Rmd `01_process_gwas_summary.Rmd` is required for processing data that comes as text files (i.e. GWAS summary stats) from the IEU GWAS pipeline or other sources. This script has to be run to convert raw data into `outcome` data frames and to extract instruments from each GWAS (in `exposure` format) and save them to be used directly in MR analysis in subsequent scripts. The names of raw files, tidy outcome data frames, and exposure instruments are all get saved in the metadata file `data_lookup.csv` upon generation. (NB the metadata file has to contain raw file names and the desired output prefixes before running this Rmd).

3. Rmd `02_mr_BMI_to_mediators.Rmd` runs univariable MR of Childhood and Adult BMI on all mediators specified in metadata file (`data_lookup.csv`). The code has to be run interactively per trait category. The results merged by trait category will be stored in `Results` directory outside the codebase. After the analysis, forest plots can be created for each trait category. To recreate the plots, don't need to rerun the full analysis, can just read in the merged files. The plots will be saved in the codebase in `figures/`. 


4. Rmd `03_mr_mediators-to-BC.Rmd` is used to run univariable MR of all mediators (and BMI) on Breast cancer (`ieu-a-1126`). The code has to be run interactively per trait category, the results are stored in `Results` outside the codebase. After the analysis, forest plots can be created for each trait category. To recreate the plots, don't need to rerun the full analysis, can just read in the merged files. The plots will be saved in the codebase in `figures/`. 

5. Rmd `04_mvmr_run_analysis.Rmd` runs four types of MVMR with each mediator, but first we run MVMR with BMIs only (Analysis 0).  

	*	 (Analysis 0) multivariable MR: childhood BMI + adult BMI -> Breast Cancer 
	*	 (Analysis 1) multivariable MR: childhood BMI + adult BMI -> mediator
	*	 (Analysis 2) multivariable MR: childhood BMI + mediator -> Breast Cancer 
	*	 (Analysis 3) multivariable MR: childhood BMI + adult BMI + mediator -> Breast Cancer 
	*	 (Analysis 4) multivariable MR: childhood BMI + childhood height + mediator -> Breast Cancer 


	The code has to be run interactively per trait category, the results are stored in `Results/trair_category/mediator/` outside the codebase. The analysis is structured as a large for-loop that will perform all four MVMR for each mediator in the selected trait category when individually specified T/F outside the loop. After all analyses have been performed, the mediators within each trait category are collated into a single dataframe and saved in `Results/trait_category/merged/` directory.

  MVMR analysis is performed using modified code from `TwoSampleMR` package and also `MVMR` package for comparison. MVMR sensitivity tests are done using `MVMR` package and phenotypic correlations estimated using `metaCCA` package (see `code_supplementary/create_pheno_covariance_matrix.Rmd`)

6. Rmd `05_mvmr_create_plots.Rmd` creates forest plots from all MVMR analyses for each trait category separately and by trait category. The plots are saved to `figures` within the codebase. The code also creates summary plot of all direct Childhood BMI estimates from Analyses 2 & 3.


7. Rmd `06_mediation_analysis.Rmd` contains a workflow for performing MR mediation analysis using Difference and Product method, with CIs calculation using Delta and Propagation of Errors methods. The script also contains methods for plotting the calculated indirect estimates (+CIs) as forest plots. 


8. **(WIP)** Rmd `07_generate_MR_reports.Rmd` calls report generation Rmd `MR_report.Rmd`:  
	
	For each mediator independently, this Rmd runs the following analyses: 
	*	 univariable MR: BMI (childhood/adult) -> mediator
	*	 univariable MR: mediator -> Breast Cancer 
	*	 multivariable MR: childhood BMI + adult BMI -> mediator
	*	 multivariable MR: childhood BMI + mediator -> Breast Cancer 
	*	 multivariable MR: childhood BMI + adult BMI + mediator -> Breast Cancer 
	*	 mediation analysis calculations (product, difference)
	*	 sensitivity analyses for univariable MR
	*	 sensitivity analyses for multivariable MR 

	The report is saved as a PDF with all sensitivity plots and tables with estimated from all MR analyses. 



_From supplementary_

`extract_from_GWAScatalog.Rmd` was used for looking up traits (breast size and mammographic density) in GWAS catalogue to extract the SNPs to use as instruments fro this traits. However, for most studies it ended up bbeing easier (more reliable) to just get the data from supplementary information provided with the papers. The code loads raw data, cleans it, and saves in the format that is ready to be used as instruments in MR analysis.

`create_pheno_covariance_matrix.Rmd` contains the workflow for applying the phenotypic correlation function from `metaCCA` package to calculate pheno_cor scores between any two GWAS traits. Tho pheno_cor score is required for MVMR sensitivity tests (F-statistics calculation). The process is time intensive, so all pheno_cor values needed for this project have been pre-calculated and are now stored in `metadata/pheno_correlations.csv`.

_From exploratory_

`experimental_MVMR.Rmd` is an adhoc testing script for other MVMR packages, supl. sensitivity analyses, as well as different combinations of exposure/mediators 3+ MVMR analyses. 























[1]: https://uob-my.sharepoint.com/:x:/r/personal/ny19205_bristol_ac_uk/Documents/Documents%20-%20OneDrive/Mini-project2/project_bibliography.xlsx?d=wff10dc7fdc9f4e0b9bc40326f7c43cc9&csf=1&web=1&e=QgyR29
[2]: https://uob-my.sharepoint.com/:x:/g/personal/ny19205_bristol_ac_uk/EX5_32z2ksVAthk1hnaUozsBfQg3QnXWEkT-5yOeXnO8SQ?e=bO4hOE
[3]: https://uob-my.sharepoint.com/:f:/g/personal/ny19205_bristol_ac_uk/EgGuNOAusytCo0f8Ods2NwABbHLq5hCN_xe5ZhyCgXBGyg?e=f0Xked

