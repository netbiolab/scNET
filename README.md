# scNET

Recontructing celltype-specific networks from single-cell RNA sequencing data
## WorkFlow introduction
![](image/introduction.png)
To construct a cell-type specific gene-gene co-expression network using single-cell RNA sequencing data, it is necessary to address the issue of sparsity. To achieve this, scNET uses bigSCale2, SAVER, and SuperCell to pre-process the count matrix, and then calculate the Pearson correlation coefficient between genes. After inferring co-expression links through 3 methods, scNET converts the PCC values into LLS and the three co-expression networks were integrated through WS score defined as the following equation: WS=〖LLS〗_0+ ∑_(i=1)^n▒〖LLS〗_i/(W×i)  ,for all LLS ≥T where 〖LLS〗_0 indicates the maximum LLS for the gene pairs and 〖LLS〗_i are sorted LLS scores by decreasing order. The weight factor, W, and LLS threshold, T, are optimized for maximizing the area under the plot of LLS versus gene coverage. The WS score for each gene pair is then re-scored using the same LLS scheme.

## Requirements
## Run with example data
