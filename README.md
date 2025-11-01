# Intestinal organoid zoo of 26 species reveals conserved and divergent programs of biology

This repository accompanies our study establishing a **“organoid zoo”** of adult stem cell–derived intestinal organoids from 26 vertebrate species, providing a unified experimental and analytical platform to explore regeneration, evolution, and host–pathogen interactions across mammals and birds.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17501594.svg)](https://doi.org/10.5281/zenodo.17501594)
---

## Overview

Biodiversity represents evolutionary strategies for survival and adaptation.  
While genomic sequencing has expanded to diverse species, cell-based comparative biology has lacked a unified platform.  
Here, we introduce a scalable framework for comparative organoid research — enabling systematic analyses across species in a standardized medium.

Our study demonstrates that this platform supports multiple applications:

- Cross-species intestinal organoid culture under standardized medium conditions
- Single-cell RNA sequencing (scRNA-seq) and comparative transcriptomics  
- Genome editing and in vitro regeneration assays  
- Viral infection modeling across mammals and avians  
- Drug screening and phenotype evaluation 

---

## Contents

This repository includes:

- Core analysis scripts and notebooks for scRNA-seq processing and analyzing 
- Code for gene expression comparison  
- Processed count matrices and metadata
- Environment lock files (`renv.lock`, `environment.yml`) for full reproducibility  

---

## Reproducibility

This repository provides both R and Python environments to reproduce all analyses.

### R Environment
```r
install.packages("renv")
renv::restore()