# Scripts for "Bioinformatic Analysis of Temporal and Spatial Proteome Alternations During Infections"
Technological advances in mass spectrometry (MS)-based proteomics and bioinformatics allow achieving a temporal and spatial resolution of the infection process at previously unseen levels. The typical output of a quantitative MS experiment that maps temporal and/or spatial changes during an infection includes highly complex, multi-dimensional data matrices with protein abundances across space or time, represented by ion intensities or spectral counts, depending on the MS approach. Such data are challenging to analyze and interpret. This repository provides machine learning and visualization scripts for the analysis of temporal and spatial proteomics data. The scripts included in this repository can also be applied to many other types of proteomic datasets. 

# The Repository Contains:
1. Scripts to analyze and generate temporal figures
- [Code:Temporal analysis](https://github.com/Babulab-bioc/TempSpac/blob/main/R/Unsupervised_Temporal.R)
2. Scripts to analyze and generate spatial figures
- [Code:Spatial analysis (unsupervised learning)](https://github.com/Babulab-bioc/TempSpac/blob/main/R/Unsupervised_Spatial.R)
- [Code:Spatial analysis (supervised learning)](https://github.com/Babulab-bioc/TempSpac/blob/main/R/ml_learning_spatial.R)

> _Note: spatial proteomics data analyses pipelines provided here use caret and factoextra R packages._


Kuhn, M., Wing, J., Weston, S., Williams, A., Keefer, C., Engelhardt, A., et al. (2020). caret: Classification and regression training. R package version 6.0–86. https://CRAN.R-project.org/package=caret [Accessed March 21, 2020].

Kassambara, A., and Mundt, F. (2020). factoextra: extract and visualize the results of multivariate data analyses. R Package version 1.0.7. https://CRAN.R-project.org/package=factoextra [Accessed April 2, 2020].

# Sample Data Description:
1. A demo spatial proteomics dataset is from **Proteome In Space and Time during Cytomegalovirus Infection. Cell Syst. 3, 361-373.** doi:10.1016/j.cels.2016.08.012. **Table S2(uninfected sample at time point 96)**.

2. A demo temporal proteomics dataset is  from the [recount project](https://jhubiostatistics.shinyapps.io/recount/) with Sequence Read Archive (SRA) accession **SRP049355**.


# License
This project is licensed under the MIT License - see the LICENSE.md file for more details.

If using these scripts in your data analyses pipelines, please cite our paper: 

**Rahmatbakhsh, M., Gagarinova, A., and Babu, M. (2021). Bioinformatic Analysis of Temporal and Spatial Proteome Alternations During Infections. Submitted to _Frontiers in Genetics, Accepted._**
