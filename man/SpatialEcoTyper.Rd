% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialEcoTyper.R
\name{SpatialEcoTyper}
\alias{SpatialEcoTyper}
\title{Identify Spatial EcoTypes from Single-cell Spatial Data (A Single Sample)}
\usage{
SpatialEcoTyper(
  normdata,
  metadata,
  outprefix = "SE",
  radius = 50,
  resolution = 0.5,
  nfeatures = 3000,
  min.cts.per.region = 1,
  npcs = 20,
  k.sn = 50,
  k = 20,
  min.cells = 5,
  min.features = 10,
  iterations = 5,
  minibatch = 5000,
  ncores = 1,
  grid.size = round(radius * 1.4),
  filter.region.by.celltypes = NULL
)
}
\arguments{
\item{normdata}{A matrix representing normalized gene expression data,
where rows correspond to genes and columns correspond to cells.}

\item{metadata}{A data frame containing metadata associated with each cell.
Must include spatial coordinates (e.g., X and Y) as well as cell type annotations.
The row names of the \code{metadata} must match the column names of the \code{normdata}.}

\item{outprefix}{Character string specifying the prefix for output file names.}

\item{radius}{Numeric specifying the radius (in the same units as spatial coordinates)
for defining spatial neighborhoods around each cell. Default is 50.}

\item{resolution}{Numeric specifying the resolution for Louvain clustering (default: 0.5).}

\item{nfeatures}{Integer specifying the number of top variable features (default: 3000) used for PCA.}

\item{min.cts.per.region}{Integer specifying the minimum number of cell types required for a microregion.}

\item{npcs}{Integer specifying the number of principal components (PCs) (default: 20).}

\item{k.sn}{Integer specifying the number of spatial nearest neighbors (default: 50) for constructing similarity network.}

\item{k}{Integer specifying the number of spatial nearest neighbors (default: 20) used to construct spatial meta-cells.}

\item{min.cells}{Minimum number of cells / spatial-meta-cells (default: 5)  expressing a feature/gene.}

\item{min.features}{Minimum number of features (default: 10) detected in a cell / spatial-meta-cell.}

\item{iterations}{Integer specifying the number of iterations (default: 5) for SNF analysis.}

\item{minibatch}{Integer specifying the number of columns to process in each minibatch in the SNF analysis.
Default is 5000. This option splits the matrix into smaller chunks (minibatch), thus reducing memory usage.}

\item{ncores}{Integer specifying the number of CPU cores to use for parallel processing.}

\item{grid.size}{Numeric specifying the grid size for spatial discretization of coordinates. By default,
this size is determined based on the specified radius (radius*1.4 µm). Increasing the grid.size will
downsample microregions and expedite the analysis, while it might eliminate cells located between bins
from the SE discovery analysis.}

\item{filter.region.by.celltypes}{A character vector specifying the cell types to include in the analysis.
Only spatial microregions that contain at least one of the specified cell types will be analyzed, while regions
lacking these cell types will be excluded from the SE discovery process. If NULL, all spatial microregions will
be included, regardless of cell type composition.

#' @return A list containing two elements:
\describe{
  \item{obj}{A seurat object constructed from fused similarity network of sptial microregions}
  \item{metadata}{Updated \code{metadata}, with a new column (`SE`) added}
}}
}
\description{
This function identifies spatially distinct cellular ecosystems (SE) from a single sample.
}
\examples{
# See https://digitalcytometry.github.io/spatialecotyper/docs/articles/SingleSample.html
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(googledrive))
suppressPackageStartupMessages(library(R.utils))
library(SpatialEcoTyper)

drive_deauth() # Disable Google sign-in requirement
drive_download(as_id("13Rc5Rsu8jbnEYYfUse-xQ7ges51LcI7n"), "HumanMelanomaPatient1_subset_counts.tsv.gz")
drive_download(as_id("12xcZNhpT-xbhcG8kX1QAdTeM9TKeFAUW"), "HumanMelanomaPatient1_subset_scmeta.tsv")

# Load single-cell gene expression matrix. Rows are genes, columns are cells.
scdata <- fread("HumanMelanomaPatient1_subset_counts.tsv.gz",
                sep = "\t",header = TRUE, data.table = FALSE)
rownames(scdata) <- scdata[, 1]  # set genes as row names
scdata <- as.matrix(scdata[, -1])
normdata <- NormalizeData(scdata)
head(normdata[, 1:5])

# Load single-cell metadata. Three columns are required, including X, Y, and CellType. Others are optional.
scmeta <- read.table("HumanMelanomaPatient1_subset_scmeta.tsv",
                     sep = "\t",header = TRUE, row.names = 1)
scmeta <- scmeta[colnames(scdata), ] # match the cell ids in scdata and scmeta
head(scmeta)

se_results <- SpatialEcoTyper(normdata, scmeta,
                              outprefix = "Melanoma1_subset",
                              radius = 50, ncores = 2)
# Extract the Seurat object and updated single-cell metadata
obj <- se_results$obj # A Seurat object
obj
scmeta <- se_results$metadata \%>\% arrange(SE) # Updated single-cell meta data, with SE annotation added
head(scmeta)
}
