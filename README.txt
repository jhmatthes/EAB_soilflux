The RStudio project directory EAB_soilflux.Rproj contains all of the code to reproduce the analysis from the raw data for the manuscript: 

Matthes, J.H., A.K. Lang, F.V. Jevon, S.J. Russell. Tree stress and mortality from emerald ash borer does not systematically alter short-term soil carbon flux in a mixed northeastern U.S. forest. Submitted to Forests, 21 Dec 2017.

The file EAB_soilflux_workflow.R steps through the analysis presented within the manuscript using functions with the R/ directory and data from the input/ directory. The workflow produces output files and figures that match those presented within the manuscript to the output/ directory.


clean_experiment_data directory:
Contains the calculated soil CO2 and CH4 flux measurements and associated covariates that will be submitted to the Dryad data repository upon publication of the manuscript. 


input directory: 
climate_NOAA: weather data from the Concord Airport (GHCND:USW00014745) station from 1981-2010 (climatology data) and 2015-2017 (experiment data)
raw_CO2CH4_data: raw data files from the Los Gatos Research Ultraportable Greenhouse Gas Analyzer and the Picarro GasScouter analyzer used for soil-atmosphere trace gas measurements
TF_insttempmoisture.csv: instantaneous soil moisture and soil temperature data collected at the time of flux measurement
TF_SoilCovariates.csv: soil covariates (pH, texture, root biomass) measured at the conclusion of the experiment
TF_TreeHealth.csv: EAB impact status of study trees within the experiment
times_key.csv: flux measurement start times that correspond to the soil-atmosphere analyzer measurements 


output directory:
data: output tables produced by the workflow code that corresponds to the Bayesian mixed models within the manuscript
figures: figures produced from the analysis that are described within the manuscript


R directory: 
functions: set of R functions used within the workflow code to format raw data and to create figures with data and model output
stan: the Stan Bayesian mixed models code used to fit the models described within the manuscript


