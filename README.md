# scNET

Reconstruction of celltype-specific networks from single-cell RNA sequencing data
## WorkFlow introduction
![](image/introduction.png)


To construct a cell-type specific gene-gene co-expression network using single-cell RNA sequencing data, it is necessary to address the issue of sparsity. To achieve this, scNET uses bigSCale, SAVER, and SuperCell to pre-process the count matrix, and then calculate the Pearson correlation coefficient between genes. After inferring co-expression links through 3 methods, scNET converts the PCC values into LLS(Log likelihood score) and the three co-expression networks were integrated through WS score euqation. The WS score for each gene pair is then re-scored using the same LLS scheme.

## Requirements
Make conda environment that R version is 4.0, and install packages using dependencies.txt file.
```bash
conda create --name scNET R==4.0
conda activate scNET
conda install --file dependencies.txt
```
Install additional packages that cannot be installed through Conda
```R
devtools::install_github("iaconogi/bigSCale2")
install.packages("igraph")
install.packages("RANN")
install.packages("WeightedCluster")
install.packages("corpcor")
install.packages("weights")
install.packages("Hmisc")
install.packages("Matrix")
install.packages("patchwork")
install.packages("plyr")
install.packages("irlba")
devtools::install_github("GfellerLab/SuperCell")
```
## Usage
The core function the run the entire pipeline is automated with scNET.sh. Briefly, it takes a seurat object or a count data file (tsv or csv), and with specific parameter outputs a inferred network in the form of edge list (gene, gene, LLS weight). prefix can be set for example as "CRC_CD8T" if the input single cell data is CD8+ T cells from colorectal cancer samples. calculated bigSCale object and SAVER-imputed values can be separately saved if BSreuse and SVreuse is set to TRUE.
```bash
scNET.sh [input] [filename] [outdir] [BScore] [SVcore] [BSreuse] [SVreuse] [package path]
```
## Parameters
[input] : full path of input gene expression matrix( it can be seurat object or gene-by-cell matrix )<br/>
[filename] : output file prefix name<br/>
[outdir] : output directory<br/>
[BScore] : number of core used for BigSCale<br/>
[SVcore] : number of core used for SAVER<br/>
[BSreuse] : Whether users reuse bigSCale matrix<br/>
[SVreuse] : Whether users reuse SAVER matrix<br/>
[package path] : directory path of this repository<br/>

We provide an examplary dataset of CD8+ Tcells from CRC context (Qian, J., Olbrecht, S., Boeckx, B. et al. https://doi.org/10.1038/s41422-020-0355-0) in the example/input folder. by setting an output folder as ./example_output/ files with CRC_T prefix will be generated within that folder. the CRC_T_scNET_Final_network.tsv is the final inferred de-novo network by combining edges inferred from each single-cell preprocessing methods and filtering them based on log-likelihood score (LLS) computed with the gold standard edge pairs provided in the input folder (GS_gold_standard_pairs.rds)

when the process finishes compare the output results with the one provided by us in example/outputfiles/...
Brief description of outputfiles if you want to add it.

```bash
scNET.sh /example/input/CRC_T_cell_count_matrix.txt CRC_T ./example_output/ 10 10 F F /scNET/
```
Setting multiple cores for BScore and SVcore parameter (10 cores each for the example code above) will reduce the running time.
## Output file
description of file | example of files|
---|---|
**preprocessed object**|CRC_T_cell_bs.rds,CRC_T_cell_sv.rds,CRC_T_cell_supercell.tsv|
**bin filtering threshold**|CRC_T_cell_bin_automation.pdf|
**network after regression**|CRC_T_cell_BS_fit,CRC_T_cell_BS_fit,CRC_T_cell_SC_fit|
**integrated network**|CRC_T_cell_final_scFGN.tsv|


## Citation
