# Allen Institute Patch-seq documents and tools
Patch-seq is powerful technique that allows for multimodal characterization of individual neurons – electrophysiological, morphological and transcriptomic. At the Allen Institute we have optimized this technique to efficiently collect high-quality data. On this GitHub repo we provide links to the manuscript describing this optimized technique and associated data, a detailed protocol, Allen Institute manuscripts that have utilized this technique, and links to Allen Institute resources and software. In addition, this repo includes an R package of [Patch-seq tools](https://github.com/AllenInstitute/patchseqtools#patch-seq-tools) for quality control and cell typing of Patch-seq cells.

## Patch-seq detailed protocol on protocols.io
https://www.protocols.io/view/patch-seq-recording-and-extraction-detailed-protoc-bw6gphbw

## Scaled, high fidelity electrophysiological, morphological, and transcriptomic cell characterization
https://elifesciences.org/articles/65482

## Multichannel Igor Electrophysiology Suite (MIES) Github repo
https://github.com/AllenInstitute/MIES

## Allen Institute Patch-seq papers
Mouse Visual Cortex https://www.sciencedirect.com/science/article/pii/S009286742031254X

Human L2/3 - https://www.nature.com/articles/s41586-021-03813-8

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

