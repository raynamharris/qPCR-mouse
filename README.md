# Integrative Hippocampus Module
This repo will contain all the data, scripts, and results from NS&B student projects form the Integrative Hippocampus Module.

## Organization

### bin
All the R scripts and other programs will be stored in bin. Additionally, the accomanying R markdown (.Rmd) and html files for project presentations will be stored here.
- SampleWrangling.R
- SampleWrangling.Rmd	
- SampleWrangling.html
- mcmc.qpcr.example.rmh


## qPCR Data Collection
qPCR data are collected using the QuantStudio3 from Thermo Fisher/ Applied Biosystems and analyzed using the QuantStudio Design and Analysis Software.
Data are collected on a USB drive then transfered to the data diretory in each scientists data directory.
The raw data are stored in .eds files that can only be read with the QuantStudio software. The data are exported as .csv files for anlaysis.

## qPCR Data Wrangling and Analysis
The .csv files are imported and wrangled in R.
The data are analyzed with the MCMC.qPCR package from the Misha Matz lab.
