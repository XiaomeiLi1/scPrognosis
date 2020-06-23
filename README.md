# scPrognosis 
## A novel single-cell based method for breast cancer prognosis
### Xiaomei Li<sup>1</sup>, Lin Liu<sup>1</sup>, Greg Goodall<sup>2,3</sup>, Andreas Schreiber<sup>2</sup>, Taosheng Xu<sup>4</sup>, Jiuyong Li<sup>1</sup>, Thuc Duy Le<sup>1</sup>

1. UniSA STEM, University of South Australia, Mawson Lakes, SA, Australia
2. Centre for Cancer Biology, an alliance of SA Pathology and University of South Australia, Adelaide, SA, Australia
3. School of Medicine, Discipline of Medicine, University of Adelaide, SA, Australia
4. Institute of Intelligent Machines, Hefei Institutes of Physical Science, Chinese Academy of Sciences, Hefei, 230031, China

This repository includes the scripts and data of the proposed method. 

The scripts are in the R folder and include:

- benchdb.R - Script for batch runing on the datasets
- benchstepwise.R - Script for stepwise Cox model (forward)
- buildNet.R - Script for building dynamic gene co-express network
- dataForWanderlust.R - Script for preparing data for the Wanderlust method
- DownstreamAnalysis.R - Script for GO enrichment analysis and KEGG pathway enrichment analysis
- geneRank.R - Script for obtaining gene rank based on the integrative model
- getPseudotime.R - Script for getting the normalized pseudotime
- grid.search.R - Script for search parameters of the integrative model
- indepTest.R - Script for independent tests
- main.R - Main script for the proposed method
- survModels.R - Script for wraped survival model
- utils.R - Script for other support functions

The input data are in the data folder and include:
- EMTmarkers.rda - EMT gene markers
- README.txt - Deascription of data
- Benchmark.signatures.rda - The gene signatures of six benchmark methods
- pseudoTime.rda - VIM-time and EMT-time used in the paper
- seed.rda - The seeds used in the experiment to repreduce the results

The output data are in the output folder and include:
- Output R files of scPrognosis and benchmark methods
- Overlap of signatures for breast cancer prognosis
- The grid-search results of scPrognosis

Notes:

(1) The single-cell data (GSE114397), GEO, TCGA753, UK can be downloaded from the link(https://github.com/XiaomeiLi1/scPrognosis/releases/tag/V1.0.1). The METABRIC data need to be download from the EMBL-EBI repository (https://www.ebi.ac.uk/ega/, accession number EGAS00000000083, require individual access agreement).

(2) The source code of Wanderlust can be downloaded from http://www.c2b2.columbia.edu/danapeerlab/html/wanderlust.html.
