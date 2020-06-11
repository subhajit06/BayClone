This is the readme file for BayClone package.
This package contains the four R code files and an example data to run the MCMC algorithm in the BayClone paper [Sengupta, 2014].

run_mcmc_simu.R:  This is the main example code to run MCMC on simulation data. It returns posterior estimates and generates heatmaps.
mcmc_cIBP.R:      This is the function of MCMC sampler for finite cIBP model.
rbetaDir.R:       This is the function to generate random beta-dirichlet variables.
compute_CPO.R:    This is the function to compute CPO for model evaluation.

simu_data.RData:  This is the example simulation data file.

In order to run the example: extract the package, go to the directory, copy paste the code in "run_mcmc_simu.R" to R console, and Z and w plots will be generated along with LPML values will be printed on the console.

Should you have any questions, please contact jinwang8@illinois.edu