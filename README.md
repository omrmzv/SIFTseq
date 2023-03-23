# SIFT-seq pipeline

**Written By:** Alexandre Cheng, Omary Mzava and Adrienne Chang from [De Vlaminck Lab](https://devlaminck.bme.cornell.edu/) at Cornell University, Ithaca, New York.



**Title:** [A metagenomic DNA sequencing assay that is robust against environmental DNA contamination](https://doi.org/10.1038/s41467-022-31654-0)
 
**Contents:** This repository contains the SIFT-seq pipeline with necessary scripts for removing environmental DNA contaminants. 
We also include R scripts required to analyze the processed datasets. 

**Implementation of SIFT-seq:**
We tag DNA by bisulfite salt-induced conversion of unmethylated cytosines to uracils. Uracils created by bisulfite treatment are converted to thymines in subsequent DNA synthesis steps that are part of DNA sequencing library preparation. After DNA sequencing, contaminating DNA introduced after tagging can then be directly identified based on the lack of cytosine conversion. We developed a bioinformatics procedure to differentiate sample-intrinsic microbial DNA, contaminant microbial DNA, and host-specific DNA after SIFT-seq tagging. Here, we present the pipeline to perform filtering.


![github](https://user-images.githubusercontent.com/62556613/173901368-d85e75fb-78f8-43b1-b8d1-5300c706442d.png)

 ## **Quick Start**
 - clone the SIFT-seq pipeline
 - Set the paths in the config/config.yaml to point to your software directory
 - Download and prepare microbe reference database
 - Running the pipeline can be done by the following line of code
 
 ```
 snakemake -j {number of cores}
 ```
 
 
