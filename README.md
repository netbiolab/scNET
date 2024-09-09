# scNET

Recontructing celltype-specific networks from single-cell RNA sequencing data
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

For example data,

```bash
scNET.sh /example/input/CRC_T_cell_count_matrix.txt CRC_T /scNET/output/ 10 10 F F /scNET/
```
Using multiple cores will reduce the running time.

## Citation
