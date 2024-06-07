# NOTE: This page is under construction. 

## Sarbecovirus disease susceptibility is conserved across viral and host models 

This document describes how to reproduce analyses in [this manuscript](https://www.biorxiv.org/content/10.1101/2023.10.11.561544v1). 

Environment prep
-----------------------

This project uses a Docker container to produce an environment similar to that used in the original analysis (e.g. R v4.2.1 and R package versions available on August 1, 2023). In order to run this container you will need [Docker](https://docs.docker.com/get-docker/) installed. 

Build the docker container:

```
docker build . -t cov 
```

Run the docker container, opening a terminal session within the container:

```
docker run -e PASSWORD=pw123 --rm -v $(pwd):/home/rstudio/work -p 8787:8787 -it cov /bin/bash
```

Navigate to the working directory: 

```
cd home/rstudio/work 
```

Prepare R/qtl files
-----------------------

Run the following code to produce a file in the format required by R/qtl.

```
mkdir -p derived_data/rqtl_files
Rscript rqtl_file_proc.R -virus SARS-CoV
Rscript rqtl_file_proc.R -virus SARS2-CoV
Rscript rqtl_file_proc.R -virus HKU3-CoV
```
These commands will generate the following .csv files, which can be imported into R and analyzed using the `Rqtl` package: 
* `derived_data/rqtl_files/SARS1_CC006xCC044_rqtl.csv` 
* `derived_data/rqtl_files/SARS2_CC006xCC044_rqtl.csv` 
* `derived_data/rqtl_files/HKU3_CC006xCC044_rqtl.csv` 


Inbred parent analysis 
-----------------------

Intercross analysis
-----------------------

Candidate gene analysis
-----------------------

