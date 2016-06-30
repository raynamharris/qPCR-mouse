# Integrative Molecular Neuroscience
## Mouse qPCR Projects

This repo contains all sorts of stuff that I am some students are working on this summer. 

We will be doing some qPCR experiments on mouse hippocampus tissues. 

We will be using the MCMC.qPCR package to analyze data.

## qPCR Data Collection
qPCR data are collected using the QuantStudio3 from Thermo Fisher/ Applied Biosystems and analyzed using the QuantStudio Design and Analysis Software.
Data are collected on a USB drive then transfered to the data diretory in each scientists data directory.
The raw data are stored in .eds files that can only be read with the QuantStudio software. The data are exported as .csv files for anlaysis.

## qPCR Data Wrangling and Analysis
The .csv files are imported and wrangled in R.
The data are analyzed with the MCMC.qPCR package from the Misha Matz lab.
