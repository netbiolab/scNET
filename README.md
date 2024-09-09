# scNET

Recontructing celltype-specific networks from single-cell RNA sequencing data
## WorkFlow introduction
![](image/introduction.png)


To construct a cell-type specific gene-gene co-expression network using single-cell RNA sequencing data, it is necessary to address the issue of sparsity. To achieve this, scNET uses bigSCale, SAVER, and SuperCell to pre-process the count matrix, and then calculate the Pearson correlation coefficient between genes. After inferring co-expression links through 3 methods, scNET converts the PCC values into LLS(Log likelihood score) and the three co-expression networks were integrated through WS score euqation. The WS score for each gene pair is then re-scored using the same LLS scheme.

## Requirements

## Usage
