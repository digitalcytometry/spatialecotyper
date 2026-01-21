#' Integrate Multiple Spatial Transcriptomics Datasets to Identify Conserved Spatial Ecotypes
#'
#' This function performs SpatialEcoTyper analysis on multiple spatial transcriptomics datasets.
#' It normalizes the input data, performs SpatialEcoTyper analysis on each dataset,
#' and integrates the results across samples.
#'
#' @param data_list A named list of expression matrices where each matrix represents
#' gene expression data for a sample. The columns of each matrix correspond to cells,
#' and the rows correspond to genes. Sample names should be used as list names.
#' Otherwise, the samples will be named as 'Sample1' through 'SampleN'.
#' @param metadata_list A named list of metadata data frames where each data frame
#' contains metadata corresponding to the cells in the expression matrices. Each row
#' should correspond to a column (cell) in the expression matrices. Each metadata
#' should include at least three columns, including X, Y and CellType.
#' @param outdir Directory where the results will be saved. Defaults to the current
#' directory with a subdirectory named "SpatialEcoTyper_results_" followed by the current date.
#' @param normalization.method Method for normalizing the expression data. Options
#' include "None" (default), "SCT", or other methods compatible with Seurat's `NormalizeData` function.
#' @param nmf_ranks Integer or a vector specifying the number of clusters (10 by default).
#' When an integer vector is supplied, the function will test all supplied numbers and
#' select the optimal number, which takes time.
#' @param nrun.per.rank An integer specifying the the number of runs per rank for NMF (default: 30).
#' @param min.coph Numeric specifying the minimum cophenetic coefficient required for a rank to be optimal.
#' @param radius Numeric specifying the radius (in the same units as spatial coordinates)
#' for defining spatial neighborhoods around each cell. Default is 50.
#' @param min.cts.per.region Integer specifying the minimum number of cell types required for a spatial neighborhood.
#' @param nfeatures An integer specifying the maximum number of top variable genes to select for each cell type.
#' @param min.features An integer specifying the minimum number of shared features (genes) required across samples.
#' @param Region Character string specifying the column name in metadata data frames
#' containing region annotations (default: NULL). Pathologist annotation is recommended if available.
#' @param downsample.by.region A Boolean indicating whether to perform downsampling
#' to balance the number of spatial clusters across regions.
#' @param subresolution Numeric specifying the resolution for clustering within each sample.
#' @param minibatch Integer specifying the number of columns to process in each minibatch in the SNF analysis.
#' Default is 5000. This option splits the matrix into smaller chunks (minibatch), thus reducing memory usage.
#' @param ncores Integer specifying the number of cores for parallel processing. Default is 1.
#' @param seed An integer used to seed the random number generator for NMF analysis.
#' @param filter.region.by.celltypes A character vector specifying the cell types to include in the analysis.
#' Only spatial neighborhoods that contain at least one of the specified cell types will be analyzed, while regions
#' lacking these cell types will be excluded from the SE discovery process. If NULL, all spatial neighborhoods will
#' be included, regardless of cell type composition.
#' @param ... Additional arguments passed to the `SpatialEcoTyper` function.
#'
#' @details This function takes a list of gene expression matrices and corresponding metadata,
#'          normalizes the data if specified, performs SpatialEcoTyper on each sample, and
#'          integrates the results across multiple samples to identify conserved spatial ecotypes.
#'
#' @return The function saves the results in the specified output directory.
#'
#' @examples
#' # See https://digitalcytometry.github.io/spatialecotyper/docs/articles/Integration.html
#'
#' @export
#'
MultiSpatialEcoTyper <- function(data_list,
                                 metadata_list,
                                 outdir = "./",
                                 normalization.method = "None",
                                 nmf_ranks = 10,
                                 nrun.per.rank = 30,
                                 min.coph = 0.95,
                                 radius = 50,
                                 min.cts.per.region = 1,
                                 nfeatures = 3000,
                                 min.features = 10,
                                 Region = NULL,
                                 downsample.by.region = TRUE,
                                 subresolution = 30,
                                 minibatch = 5000,
                                 ncores = 1,
                                 seed = 1,
                                 filter.region.by.celltypes = NULL,
                                 ...){

  outdir <- file.path(outdir, paste0("SpatialEcoTyper_results_", Sys.Date()))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if(is.null(names(data_list)) & is.null(names(metadata_list))){
    sample_names <- paste0("Sample", 1:length(data_list))
    names(data_list) <- sample_names
    names(metadata_list) <- sample_names
  }else{
    if(is.null(names(data_list))) names(data_list) = names(metadata_list)
    if(is.null(names(metadata_list))) names(metadata_list) = names(data_list)
    sample_names <- names(data_list)
  }
  idx <- unlist(lapply(sample_names, function(ss){
    ncol(data_list[[ss]]) == nrow(metadata_list[[ss]])
  }))
  if(sum(!idx)>0){
    stop("The number of cells in expression data and meta data do not match for: ",
         paste0(sample_names[!idx], collapse = ", "))
  }

  seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
  if(normalization.method!="None"){
    message(Sys.time(), " Normalize the expression matrices")
    if(normalization.method=="SCT"){
      data_list <- mclapply(data_list, function(x){
        tmpobj = CreateSeuratObject(x) %>% SCTransform(clip.range = c(-10, 10), verbose = FALSE)
        if(seurat_version<5){
          GetAssayData(tmpobj, "data")
        }else{
          tmpobj[["SCT"]]$data
        }
      }, mc.cores = ncores)
    }else{
      data_list <- lapply(data_list, NormalizeData,
                          normalization.method = normalization.method,
                          verbose = FALSE)
    }
  }

  message(Sys.time(), " Perform SpatialEcoTyper analysis for each sample")
  SpatialEcoTyper_list <- list()
  for(ss in sample_names){
    message(Sys.time(), " SpatialEcoTyper analysis for ", ss)
    SpatialEcoTyper_list[[ss]] = SpatialEcoTyper(normdata = data_list[[ss]],
                          metadata = metadata_list[[ss]],
                          outprefix = paste0(outdir, "/", ss),
                          radius = radius,
                          resolution = subresolution,
                          min.cts.per.region = min.cts.per.region,
                          minibatch = minibatch,
                          ncores = ncores,
                          filter.region.by.celltypes = filter.region.by.celltypes,
                          ...)
  }

  IntegrateSpatialEcoTyper(SpatialEcoTyper_list = SpatialEcoTyper_list,
                           data_list = data_list,
                           outdir = outdir,
                           normalization.method = "None",
                           nmf_ranks = nmf_ranks,
                           nrun.per.rank = nrun.per.rank,
                           min.coph = min.coph,
                           nfeatures = nfeatures,
                           min.features = min.features,
                           Region = Region,
                           downsample.by.region = downsample.by.region,
                           subresolution = subresolution,
                           minibatch = minibatch,
                           ncores = ncores,
                           seed = seed)
}

