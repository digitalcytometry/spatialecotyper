
<p align="center">
<img width="800" src="logo.png">
</p>

# Deciphering spatial ecotypes from spatial, single-cell, and bulk transcriptomic data

**Spatial EcoTyper** is a versatile framework for identifying spatially distinct multicellular communities, termed spatial ecotypes, from single-cell spatial transcriptomics data. In addition, it provides unified methods for the recovery of spatial ecotypes across multiple data modalities, including spatial transcriptomics, single-cell RNA-seq, and bulk transcriptomic datasets.

## Overview of Spatial EcoTyper

**Spatial EcoTyper** is available as an R package, with comprehensive documentation accessible at https://digitalcytometry.github.io/spatialecotyper.

We provide eight comprehensive tutorials illustrating the key functionalities of Spatial EcoTyper framework:

-   **Tutorial 1**: [Discovering Spatial Ecotypes from a Single Spatial Transcriptomics Sample](https://digitalcytometry.github.io/spatialecotyper/articles/SingleSample.html)
-   **Tutorial 2**: [Discovering Conserved Spatial Ecotypes Across Multiple Spatial Transcriptomics Samples](https://digitalcytometry.github.io/spatialecotyper/articles/Integration.html)
-   **Tutorial 3**: [Identifying SE-Specific Cell States via Leave-One-Sample-Out Cross-Validation](https://digitalcytometry.github.io/spatialecotyper/articles/Discovery_SE_CellStates.html)
-   **Tutorial 4**: [NMF Model Development for Spatial Ecotype Recovery from Single-Cell and Spatial Transcriptomics Data](https://digitalcytometry.github.io/spatialecotyper/articles/TrainRecoveryModel.html)
-   **Tutorial 5**: [Recovering Spatial Ecotypes from Single-Cell Spatial Transcriptomics Data](https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_scST.html)
-   **Tutorial 6**: [Recovering Spatial Ecotypes from Single-Cell RNA-seq Data](https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_scRNA.html)
-   **Tutorial 7**: [NMF Model Development for Spatial Ecotype Deconvolution from Bulk Gene Expression Data](https://digitalcytometry.github.io/spatialecotyper/articles/TrainDeconvModel.html)
-   **Tutorial 8**: [Inferring Spatial Ecotype Abundances from Bulk Gene Expression Data](https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_Bulk.html)

**Note**: Spatial EcoTyper depends extensively on Seurat for key processes like dimensionality reduction, UMAP embedding, clustering, and visualization. Initially developed using Seurat v4.3, the tool has been thoroughly tested and validated with Seurat v5. Although UMAP embeddings and clustering results show slight differences between Seurat v4 and v5, the overall consistency remains strong, ensuring that core biological insights are preserved across both versions.

## System requirements

- **R**: Version 4.0 or higher is required.
- **Operating Systems**: This package is compatible with all operating systems and has been tested on the following platforms:
  - macOS: Big Sur, Monterey, Ventura, Sonoma, Sequoia (15.2)
  - Linux: CentOS 7.2 and High-Performance Computing (HPC) clusters

## Installation

**Spatial EcoTyper** is available as an R package and can be installed via the `BiocManager` package directly from the R console.

``` r
if(!"BiocManager" %in% installed.packages()){
  install.packages("BiocManager")
}

## Install dependencies
BiocManager::install(c("remotes", "Seurat", "NMF", "dplyr", "tidyr", "pals",
                       "parallel", "data.table", "ComplexHeatmap", 
                       "glmGamPoi", "immunogenomics/presto"))

## Install SpatialEcoTyper
BiocManager::install("digitalcytometry/spatialecotyper")
```


<details><summary>Troubleshooting dependency installation</summary>

* ERROR: dependency ‘GetoptLong’ is not available for package ‘ComplexHeatmap’

  If the installation within R console fails, you can try installing the necessary packages via `conda install` or `mamba install`.
  
  ```bash
  conda install bioconda::bioconductor-complexheatmap
  ```


* Failed to install 'presto' from GitHub: HTTP error 401. Bad credentials

  To resolve this issue, you’ll need to authenticate using a personal access token (PAT). You can generate a GitHub personal access token following the [GitHub's documentation](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens). After generating the token, set it as an environment variable in R using the following code. Replace "YOUR_TOKEN" with your actual token:
  ```r
  ## Set the token in your R environment:
  Sys.setenv(GITHUB_PAT="YOUR_TOKEN")
  
  ## Install the package from GitHub:
  BiocManager::install("immunogenomics/presto")
  ```

</details>


<details><summary>Install SpatialEcoTyper from source code</summary>
The source code of **Spatial EcoTyper** is available at https://github.com/digitalcytometry/spatialecotyper. After downloading the package, you can install it from the source code using the command:

``` r
install.packages("SpatialEcoTyper.tar.gz", repos = NULL)
```

</details>


## Contribution

If you encounter any bugs or have suggestions for improvements, please feel free to open an [issue](https://github.com/digitalcytometry/spatialecotyper/issues) or submit a pull request. Your feedback and contributions help us make the tool better for everyone.

## License
Please see the <a href="LICENSE.html" target="_blank">LICENSE</a> file.


## Citation
If you use Spatial EcoTyper, please cite: 

Wubing Zhang<sup>\*</sup>, Erin L. Brown<sup>\*</sup>, Abul Usmani<sup>\*</sup>, Noah Earland, Minji Kang, Chibuzor Olelewe, Anushka Viswanathan, Pradeep S. Chauhan, Chloé B. Steen, Hyun Soo Jeon, Susanna Avagyan, Irfan Alahi, Nicholas P. Semenkovich, Janella C. Schwab, Chloe M. Sachs, Faridi Qaium, Peter K. Harris, Qingyuan Cai, Andrew J. Gentles, James Knight, Rondell P. Graham, Antonietta Bacchiocchi, Peter C. Lucas, Ryan C. Fields, Mario Sznol, Ruth Halaban, David Y. Chen, Aadel A. Chaudhuri<sup>†</sup> and Aaron M. Newman<sup>†</sup>. **Non-invasive profiling of the tumour microenvironment with spatial ecotypes.** *Nature*, 2026.  [doi.org/10.1038/s41586-026-10452-4](https://www.nature.com/articles/s41586-026-10452-4).
