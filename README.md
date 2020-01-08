CIM-seq publication
================

### Introduction

Facilitates reproduction of figures associated with the CIM-seq publication.

### Prerequisites

Install R, git, and docker first if not installed.

### Instructions

1.  Clone the repo<br/> `git clone https://github.com/EngeLab/CIMseq.publication.git`

2.  Change to the repo directory<br/> `cd CIMseq.publication`

3.  Get the docker image<br/> `docker pull engelab/cim-seq-publication`

4.  Generate the processed data files<br/> `docker run -v $PWD:/home/CIMseq.publication engelab/cim-seq-publication:latest Rscript -e "future::plan(future::multiprocess); source('./inst/analysis/runAnalysis.R')"`

Figures are located in inst/analysis/name\_of\_the\_analysis/figures.<br/> Counts data for the sorted multiplet dataset and the mouse gut dataset are in the data folder.

*Don't forget to clean up the Docker image/container/etc.*
