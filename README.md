# Allen Institute Patch-seq documents and tool

## patchseqtools
  
R functions for gene selection and analysis of Patch-seq data.  
  
`patchseqtools` includes several functions that are used for quality control and cell typing of Patch-seq cells at the Allen Institute.  Many of these functions are wrappers for functions from other libraries (see below) and we enourage proper citation of relevant tools when using those functions.  

**Specific topics include:**  
1. Assigning quality scores to each cell (mostly wrapper functions for https://github.com/PavlidisLab/patchSeqQC)  
2. Cell type clustering using CCA (mostly wrapper functions for https://satijalab.org/seurat/Seurat_AlignmentTutorial.html)  
  
 
## Installation

```
devtools::install_github("AllenInstitute/patchseqtools")
```

## Library use cases

No vignettes currently available.  

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/patchseqtools/blob/master/LICENSE

## Level of Support

We are not planning to update the R library component of this repo unless bugs are found.  We may make occasional updates to other components of the repository with no fixed schedule.  However, we welcome community input through both issues and pull requests.

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/patchseqtools/blob/master/CONTRIBUTION

