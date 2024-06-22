## Sarbecovirus disease susceptibility is conserved across viral and host models 

This document describes how to reproduce analyses in this manuscript: 10.1016/j.virusres.2024.199399

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
Rscript rqtl_file_proc.R --args --virus=SARS-CoV
Rscript rqtl_file_proc.R --args --virus=SARS2-CoV
Rscript rqtl_file_proc.R --args --virus=HKU3-CoV
```
These commands will generate the following .csv files, which can be imported into R and analyzed using the `Rqtl` package: 
* `derived_data/rqtl_files/SARS1_CC006xCC044_rqtl.csv` 
* `derived_data/rqtl_files/SARS2_CC006xCC044_rqtl.csv` 
* `derived_data/rqtl_files/HKU3_CC006xCC044_rqtl.csv` 


Inbred parent analysis 
-----------------------

Run the following code to analyze data from the initial and repeated Collaborative Cross screens and generate plots that comprise Figure 1. 

```
Rscript parent-analysis.R
```

Several plots will be produced in `figures/screen/`.

Intercross analysis
-----------------------

Run the following code to perform QTL mapping in each infection group and produce a PDF of results. 

```
Rscript -e 'library(rmarkdown); rmarkdown::render("vr-Rqtl_SARS1.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("vr-Rqtl_SARS2.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("vr-Rqtl_HKU3.Rmd", "html_document")'
```

Candidate gene analysis
-----------------------

Run the following command to decompress the *HrS43* variant data:

```
unzip source_data/chr9-variants.zip -d source_data/
```

Run the following Rscript to perform candidate gene analysis on *HrS43*. Since this region was identified in a similar cross between CC011 and CC074 [(Schäfer et al, 2022)](https://journals.asm.org/doi/10.1128/mbio.01454-22), this code looks for genes with variants segregating between the strains in both crosses, as well as genes with variants segregating between the strains in either cross. It produces two files: `results/candidate-genes.csv` contains all candidate genes, and `results/filtered-cand-genes.csv` contains only candidate genes that are segregating between both crosses. 

```
Rscript HrS43_cand_gene_analysis.R
```

Please reach out to elrisemberg@gmail.com with any questions or concerns. 