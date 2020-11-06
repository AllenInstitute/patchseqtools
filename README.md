# Allen Institute Patch-seq documents and tools
## Patch-seq detailed protocol on protocols.io
https://doi.org/10.17504/protocols.io.bpbuminw

## High fidelity electrophysiological, morphological, and transcriptomic cell characterization using a refined Patch-seq protocol
https://doi.org/10.1101/2020.11.04.369082

## Metadata for manuscript##

## Multichannel Igor Electrophysiology Suite (MIES) Github repo
https://github.com/AllenInstitute/MIES

## Allen Institute Patch-seq papers
Mouse Visual Cortex https://doi.org/10.1101/2020.02.03.932244

Human L2/3 - https://doi.org/10.1101/2020.03.31.018820

## Allen Institute cell types web products
Morpho-Electric database http://celltypes.brain-map.org/

Single-cell RNA-seq database https://portal.brain-map.org/atlases-and-data/rnaseq

Patch-Seq database *Coming soon*

## Patch-seq tools
  
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

