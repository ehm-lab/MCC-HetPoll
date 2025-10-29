# Differential mortality risks associated to PM2.5 components: a multi-country multi-city study.

[![DOI](https://zenodo.org/badge/262352027.svg)](https://zenodo.org/doi/10.5281/zenodo.11210101)

R code and results attached to the publication:

Masselot P, et al. (2022) Differential mortality risks associated to PM2.5 components: a multi-country multi-city study. *Epidemiology*. **33-2**. [DOI:10.1097/EDE.0000000000001455](https://doi.org/10.1097/EDE.0000000000001455)

### Data

Mortality data are not available currently due to restricted data sharing agreement between the collaborators of this study. Therefore, **only scripts 3 to 5 are fully reproducible**.

The *Data* folder contains data needed to run scripts 3 to 5, *i.e.* results from the first-stage city-specific models and PM composition data.

### R code

The R code to reproduce the analysis is available. Scripts are meant to be executed in order:

- *0_PrepData.R* Loads and and links mortality, PM2.5, composition and city-specific characteristics. Performs city selection.
- *1_DataSummary.R* Creates Table 1 providing descriptive statistics.
- *2_FirstStage.R* Runs the first-stage model on each selected city to produce city-level relative risks.
- *3_SecondStage_MetaRegComposition.R* Pools first-stage results in a meta analysis including city-specific compositions and socio-economic indicators. 
- *4_Plots.R* Produces the figures of the paper.
- *5_SupplementaryResults.R* Produces supplementary results and eFigures.

### Results

Figures and Tables generated for the publication (including supplemental ones) are available in the Results folder of this repository. 

### Acknowledgements

This work was supported by the Medical Research Council of UK (Grant ID: MR/M022625/1), the Natural Environment Research Council of UK (Grant ID: NE/R009384/1), and the European Union’s Horizon 2020 Project Exhaustion (Grant ID: 820655).
