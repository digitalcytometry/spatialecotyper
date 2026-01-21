#' Identify Spatial EcoTypes from Single-cell Spatial Data (A Single Sample)
#'
#' This function identifies spatially distinct cellular ecosystems (SE) from a single sample.
#'
#' @param normdata A matrix representing normalized gene expression data,
#' where rows correspond to genes and columns correspond to cells.
#' @param metadata A data frame containing metadata associated with each cell.
#' Must include spatial coordinates (e.g., X and Y) as well as cell type annotations.
#' The row names of the \code{metadata} must match the column names of the \code{normdata}.
#' @param outprefix Character string specifying the prefix for output file names.
#' @param radius Numeric specifying the radius (in the same units as spatial coordinates)
#' for defining spatial neighborhoods around each cell. Default is 50.
#' @param resolution Numeric specifying the resolution for Louvain clustering (default: 0.5).
#' @param nfeatures Integer specifying the number of top variable features (default: 500) used for the analysis.
#' @param min.cts.per.region Integer specifying the minimum number of cell types required for a spatial neighborhood.
#' @param npcs Integer specifying the number of principal components (PCs) (default: 20).
#' @param min.cells Minimum number of cells / spatial-meta-cells (default: 5)  expressing a feature/gene.
#' @param min.features Minimum number of features (default: 10) detected in a cell / spatial-meta-cell.
#' @param iterations Integer specifying the number of iterations (default: 10) for SNF analysis.
#' @param minibatch Integer specifying the number of columns to process in each minibatch in the SNF analysis.
#' Default is 5000. This option splits the matrix into smaller chunks (minibatch), thus reducing memory usage.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#' @param grid.size Numeric specifying the grid size for spatial discretization of coordinates. By default,
#' this size is determined based on the specified radius round(radius*1.4 µm). Increasing the grid.size will
#' downsample spatial neighborhoods and expedite the analysis, while it might eliminate cells located between bins
#' from the SE discovery analysis.
#' @param filter.region.by.celltypes A character vector specifying the cell types to include in the analysis.
#' Only spatial neighborhoods that contain at least one of the specified cell types will be analyzed, while regions
#' lacking these cell types will be excluded from the SE discovery process. If NULL, all spatial neighborhoods will
#' be included, regardless of cell type composition.
#' @param k Integer specifying the number of spatial nearest neighbors (default: 20) used to construct spatial meta-cells.
#' @param k.sn Integer specifying the number of nearest neighbors (default: 50) for constructing similarity network.
#' @param dropcell Logical. If TRUE, cells that cannot be assigned to any spatial
#' ecotype (outside the radius) will be removed from the returned metadata. Default is TRUE.
#'
#' #' @return A list containing two elements:
#' \describe{
#'   \item{obj}{A seurat object constructed from fused similarity network of sptial neighborhoods}
#'   \item{metadata}{Updated \code{metadata}, with a new column (`SE`) added}
#' }
#'
#' @import Seurat
#' @importFrom parallel detectCores
#'
#' @examples
#' # See https://digitalcytometry.github.io/spatialecotyper/docs/articles/SingleSample.html
#' suppressPackageStartupMessages(library(dplyr))
#' suppressPackageStartupMessages(library(ggplot2))
#' suppressPackageStartupMessages(library(parallel))
#' suppressPackageStartupMessages(library(Seurat))
#' suppressPackageStartupMessages(library(data.table))
#' suppressPackageStartupMessages(library(googledrive))
#' suppressPackageStartupMessages(library(R.utils))
#' library(SpatialEcoTyper)
#'
#' drive_deauth() # Disable Google sign-in requirement
#' drive_download(as_id("13Rc5Rsu8jbnEYYfUse-xQ7ges51LcI7n"), "HumanMelanomaPatient1_subset_counts.tsv.gz")
#' drive_download(as_id("12xcZNhpT-xbhcG8kX1QAdTeM9TKeFAUW"), "HumanMelanomaPatient1_subset_scmeta.tsv")
#'
#' # Load single-cell gene expression matrix. Rows are genes, columns are cells.
#' scdata <- fread("HumanMelanomaPatient1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]  # set genes as row names
#' scdata <- as.matrix(scdata[, -1])
#' normdata <- NormalizeData(scdata)
#' head(normdata[, 1:5])
#'
#' # Load single-cell metadata. Three columns are required, including X, Y, and CellType. Others are optional.
#' scmeta <- read.table("HumanMelanomaPatient1_subset_scmeta.tsv",
#'                      sep = "\t",header = TRUE, row.names = 1)
#' scmeta <- scmeta[colnames(scdata), ] # match the cell ids in scdata and scmeta
#' head(scmeta)
#'
#' se_results <- SpatialEcoTyper(normdata, scmeta,
#'                               outprefix = "Melanoma1_subset",
#'                               radius = 50, ncores = 2)
#' # Extract the Seurat object and updated single-cell metadata
#' obj <- se_results$obj # A Seurat object
#' obj
#' scmeta <- se_results$metadata %>% arrange(SE) # Updated single-cell meta data, with SE annotation added
#' head(scmeta)
#' @export
#'
SpatialEcoTyper <- function(normdata, metadata,
                            outprefix = "SE",
                            radius = 50,
                            resolution = 0.5,
                            nfeatures = 500,
                            min.cts.per.region = 2,
                            npcs = 20,
                            min.cells = 5,
                            min.features = 10,
                            iterations = 10,
                            minibatch = 5000,
                            ncores = 1,
                            grid.size = round(radius*1.4),
                            filter.region.by.celltypes = NULL,
                            k = 20,
                            k.sn = 50,
                            dropcell = TRUE){

  if(ncol(normdata)!=nrow(metadata)){
    stop("The number of cells in expression data and meta data do not match.")
  }
  tmp <- PreprocessST(normdata, metadata = metadata,
                      min.cells = min.cells,
                      min.features = min.features)
  normdata <- tmp$expdat
  metadata <- tmp$metadata

  if(!"CellType" %in% colnames(metadata)) stop("The meta data should include a column (CellType) for cell type annotations")
  if(!"X" %in% colnames(metadata)) stop("The meta data should include spatial coordinates (X and Y) of cells")
  if(!"Y" %in% colnames(metadata)) stop("The meta data should include spatial coordinates (X and Y) of cells")

  ncmeta = metadata
  ncmeta$SpotID <- paste0("X", round(ncmeta$X / grid.size), "_Y", round(ncmeta$Y / grid.size))
  ncmeta <- ncmeta %>% group_by(SpotID) %>%
    mutate(Spot.X = median(X), Spot.Y = median(Y)) %>% as.data.frame
  rownames(ncmeta) = rownames(metadata)

  ### Get spatial meta cells for each cell type
  message(Sys.time(), " Construct spatial meta cells for each cell type")
  if(max(normdata)>50) normdata <- log1p(normdata+1)
  ncem <- GetSpatialMetacells(normdata, ncmeta, k = k, radius = radius,
                              ncores = ncores)
  ### Only interested in spatial regions include a certain cell type?
  if(!is.null(filter.region.by.celltypes)){
    spots = unlist(lapply(filter.region.by.celltypes, function(x){
      idx = gsub(".*\\.+", "", colnames(ncem)) %in% filter.region.by.celltypes
      spots = gsub("\\.+.*", "", colnames(ncem)[idx])
      spots
    }))
    spots = table(spots)
    spots = names(spots)[spots==length(filter.region.by.celltypes)]
    idx = gsub("\\.+.*", "", colnames(ncem))%in%spots
    if(sum(idx)<20) stop("Fewer than 20 spatial neighborhoods contain ",
                         paste0(filter.region.by.celltypes, collapse = ", "))
    ncem = ncem[, idx]
  }
  ## spatial neighborhood level meta data
  ncmeta <- ncmeta[, colnames(ncmeta)!="CID"]
  ## If categorical, take the most frequent category, if numeric, take median.
  tmpmedian <- ncmeta %>% group_by(SpotID, Spot.X, Spot.Y) %>% summarise(across(where(is.numeric), median))
  tmp_freq <- ncmeta %>% group_by(SpotID, Spot.X, Spot.Y) %>%
    summarise(across(where(is.character) | where(is.factor) | where(is.logical), mostFrequent))
  ncmeta <- merge(tmpmedian, tmp_freq, by = c("SpotID", "Spot.X", "Spot.Y"), all = TRUE) %>% as.data.frame
  rownames(ncmeta) <- ncmeta$SpotID

  spots = unique(gsub("\\.+.*", "", colnames(ncem)))
  message("\t\tTotal spatial neighborhoods: ", length(spots),
          "\n\t\tTotal spatial meta cells: ", ncol(ncem))
  rm(normdata)

  ### Get PCs for all cell types
  message(Sys.time(), " PCA for each cell type")
  emb_list <- GetPCList(ncem, nfeatures = nfeatures, min.cells = min.cells,
                        min.features = min.features, ncores = ncores)
  rm(ncem)

  ### SNF
  message(Sys.time(), " Construct cell-type-specific similarity network")
  snlist <- GetSNList(emb_list, npcs = npcs,
                      min.cts.per.region = min.cts.per.region,
                      k = k.sn, ncores = ncores)
  rm(emb_list)
  # if(!is.null(outprefix)) saveRDS(snlist, gsub("SNF_obj", "snlist", snf_outf))

  message(Sys.time(), " Similarity network fusion")
  obj = SNF2(snlist, t = iterations, minibatch = minibatch,
             ncores = ncores, verbose = TRUE)
  rm(snlist)
  obj = rankSparse(obj)
  rownames(obj) = gsub("_", "-", rownames(obj))

  message(Sys.time(), " Create Seurat object and perform clustering analysis")
  obj <- CreateSeuratObject(obj, meta.data = ncmeta[colnames(obj), ])
  seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
  if(seurat_version>=5){
    obj[["RNA"]]$data = obj[["RNA"]]$counts
  }
  VariableFeatures(obj) = rownames(obj)
  obj <- obj %>% ScaleData(do.scale = FALSE, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:10, verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE)
  obj$SE = paste0("SE", obj$seurat_clusters)
  obj@project.name = paste0(outprefix, "_radius", radius, "_nfeatures", nfeatures, "_npcs", npcs,
                            "_k.sn", k.sn, "_k", k, "_min.cells", min.cells,
                            "_min.features", min.features,
                            "_min.cts.per.region", min.cts.per.region,
                            "_iterations", iterations, "_grid.size", grid.size)
  metadata = AnnotateCells(metadata, obj, dropcell = dropcell)
  results <- list(obj = obj, metadata = metadata)
  if(!is.null(outprefix)){
    snf_outf <- paste0(outprefix, "_SpatialEcoTyper_results.rds")
    saveRDS(results, snf_outf)
    message(Sys.time(), " The Seurat object is saved into ", snf_outf)
  }
  return(results)
}
