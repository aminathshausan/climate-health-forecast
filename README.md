# climate-health-forecast
This repository makes available the source code for the manuscript: "Climate-informed short-term forecasting of influenza incidence in Australia". The computation concerns model fits,  forecast validations and plotting of 4  Bayesian spatiotemporal models. The computation is performed using Integrated Laplace Approximation (INLA) method. All model fits are implemented using the R software and data cleaning and some plotting was conducted using python software.


# Prerequisites

To run the R codes  codes you will require R(>=3.30) and  latest version of RStudio. To run python codes you will need Jupyter Notebook and latest python version

# Directory structure
```
├── project
│   ├── data
│   ├── figs
|   ├── results
│   ├── src
│   │   ├── fit_models_noCovt.R
│   │   ├── plots.R
│   │   ├── posterior_linr_combs.R
```
- The `data` folder contains all processed data required to run the codes in the `src` folder.
- The `figs` folder is a place holder for all figures generated from various codes  in the `src` folder. 
- The `results` folder contains forecast validation values for all model fits and the optimal model.
- The `src` folder contains all source codes required to produce the results and images in the manuscript.

# How to run scripts

- Run the scripts in the `src` folder to  fit the models, generate plots.
 
