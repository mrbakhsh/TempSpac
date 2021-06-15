# Scripts for "Bioinformatic Analysis of Temporal and Spatial Proteome Alternations During Infections"
Technological advances in mass spectrometry (MS)-based proteomics and bioinformatics allow achieving a temporal and spatial resolution of the infection process at previously unseen levels. The typical output of a quantitative MS experiment that maps temporal and/or spatial changes during an infection includes highly complex, multi-dimensional data matrices with protein abundances across space or time, represented by ion intensities or spectral counts, depending on the MS approach. Such data are challenging to analyze and interpret. This repository provides machine learning and visualization scripts for the analysis of temporal and spatial proteomics data. The scripts included in this repository can also be applied to many other types of proteomic datasets. 

# The Repository Contains:
1. Scripts to analyze and generate temporal figures
- [Code:Temporal analysis](https://github.com/Babulab-bioc/TempSpac/blob/main/R/Unsupervised_Temporal.R)
2. Scripts to analyze and generate spatial figures
- [Code:Spatial analysis (unsupervised learning)](https://github.com/Babulab-bioc/TempSpac/blob/main/R/Unsupervised_Spatial.R)
- [Code:Spatial analysis (supervised learning)](https://github.com/Babulab-bioc/TempSpac/blob/main/R/ml_learning_spatial.R)



# Sample Data Description:
1. A demo spatial proteomics dataset is from Proteome In Space and Time during Cytomegalovirus Infection. Cell Syst. 3, 361-373. doi:10.1016/j.cels.2016.08.012. Table S2 (uninfected sample at time point 96).

2. A demo temporal proteomics dataset is  from the [recount project](https://jhubiostatistics.shinyapps.io/recount/) with Sequence Read Archive (SRA) accession SRP049355.


# License
This project is licensed under the MIT License - see the LICENSE.md file for more details.

If using these scripts in your data analyses pipelines, please cite our paper: 

**Rahmatbakhsh M, Gagarinova A and Babu M (2021) Bioinformatic Analysis of Temporal and Spatial Proteome Alternations During Infections. _Front. Genet_. 12:667936. doi: 10.3389/fgene.2021.667936**.
