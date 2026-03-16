Within this repository are all the files needed to replicate the analysis. 

Due to file sizes, the shapefiles in the "raw" folder have been compressed. To use the contents of these folders, unzip them to the "raw" folder. Note due to GitHub file size constraints, the boundaries_p_2021_v3 shapefile was not included in the repository. This shapefile can be found at https://www.cec.org/north-american-environmental-atlas/political-boundaries-2021/

To replicate the analysis begin by running the "ProcessLakes.R" script.  Next, run the "MCMC.script.R" to run the Markov Chain Monte Carlo simulation. The remaining R scripts may be run in any order. 
