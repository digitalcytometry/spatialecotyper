  
  
<p align="center">
<img width="800" src="logo.png">
</p>

# Deciphering spatial ecotypes from spatial, single-cell, and bulk transcriptomic data

**Spatial EcoTyper** is a versatile framework designed for the systematic identification of spatial cellular communities, termed spatial ecotypes, from single-cell spatial transcriptomics data. In addition, it provides unified methods for the recovery of spatial ecotypes across multiple data modalities, including spatial transcriptomics, single-cell RNA-seq, and bulk transcriptomic datasets.

## Overview of Spatial EcoTyper

**Spatial EcoTyper** is available as an R package, with comprehensive documentation accessible at https://digitalcytometry.github.io/spatialecotyper.

We provide five comprehensive tutorials illustrating the functionalities included in the SpatialEcoTyper R package. The first tutorial demonstrates how to identify spatial ecotypes from a single-cell spatial transcriptomics data. The second demonstrates how to identified conserved spatial ecotypes across multiple samples. The third demonstrates how to develop NMF models for recovery of spatial ecotypes from unseen data. The remaining tutorials introduce how to recover spatial ecotypes from single-cell (scRNA-seq or single-cell spatial transcriptomics) and bulk (bulk RNA-seq or Visium) gene expression profiles.

-   **Tutorial 1:** [Discovery of Spatial Ecotypes from A Single Sample](https://digitalcytometry.github.io/spatialecotyper/articles/SingleSample.html)
-   **Tutorial 2:** [Discovery of Spatial Ecotypes from Multiple Samples](https://digitalcytometry.github.io/spatialecotyper/articles/Integration.html)
-   **Tutorial 3:** [Development of NMF Models for Spatial Ecotype Recovery](https://digitalcytometry.github.io/spatialecotyper/articles/TrainRecoveryModels.html)
-   **Tutorial 4:** [Recovery of Spatial Ecotypes from Single-Cell Gene Expression Data](https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_scRNA.html)
-   **Tutorial 5:** [Recovery of Spatial Ecotypes from Bulk Gene Expression Data](https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_Bulk.html)

Data used in the tutorials can be obtained from https://drive.google.com/open?id=1En0pgY6_3_u8XK2hTSTbrA9Iouhvc4TD&usp=drive_fs.

**Note**: __Spatial EcoTyper__ depends extensively on Seurat for key processes like dimensionality reduction, UMAP embedding, clustering, and visualization. Initially developed using Seurat v4.3, the tool has been thoroughly tested and validated with Seurat v5. Although UMAP embeddings and clustering results show slight differences between Seurat v4 and v5, the overall consistency remains strong, ensuring that core biological insights are preserved across both versions.

## System requirements

This package is compatible with all operating systems and has been tested on the following platforms:

macOS: Big Sur, Monterey, Ventura, Sonoma, Sequoia (15.2)
Linux: CentOS 7.2 and High-Performance Computing (HPC) clusters


## Installation

**Spatial EcoTyper** is available as an R package and can be installed via the `BiocManager` package directly from the R console.

``` r
if(!"BiocManager" %in% installed.packages()){
  install.packages("BiocManager")
}

## Install dependencies
BiocManager::install(c("remotes", "Seurat", "NMF", "dplyr", "pals",
                       "data.table", "ComplexHeatmap", "googledrive", 
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

## Authors

Spatial EcoTyper was developed in the <a href="https://anlab.stanford.edu/" target="_blank">Newman Lab</a> by Wubing Zhang.

## Citation
If you use Spatial EcoTyper, please cite:

