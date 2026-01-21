#' Integrate Multiple Spatial Transcriptomics Datasets to Identify Conserved Spatial Ecotypes
#'
#' @param SpatialEcoTyper_list A named list of SpatialEcoTyper results, each item
#' represents a list returned from the \code{SpatialEcoTyper} function.
#' @param data_list A named list of expression matrices where each matrix represents
#' gene expression data used for the SpatialEcoTyper analysis. The list name should match
#' that of \code{SpatialEcoTyper_list}.
#' @param outdir Directory where the results will be saved. Defaults to the current
#' directory with a subdirectory named "SpatialEcoTyper_results_" followed by the current date.
#' @param normalization.method Method for normalizing the expression data. Options
#' include "None" (default), "SCT", or other methods compatible with Seurat's `NormalizeData` function.
#' @param nmf_ranks Integer or a vector specifying the number of clusters (10 by default).
#' When an integer vector is supplied, the function will test all supplied numbers and
#' select the optimal number, which takes time.
#' @param nrun.per.rank An integer specifying the the number of runs per rank for NMF (default: 30).
#' @param min.coph Numeric specifying the minimum cophenetic coefficient required for a rank to be optimal.
#' @param nfeatures An integer specifying the maximum number of top variable genes to select for each cell type.
#' @param min.features An integer specifying the minimum number of shared features (genes) required across samples.
#' @param Region A character string specifying the column name in metadata data frames
#' containing region annotations (default: NULL). Pathologist annotation is recommended if available.
#' @param downsample.by.region A Boolean indicating whether to perform downsampling
#' to balance the number of spatial clusters across regions.
#' @param subresolution A numeric specifying the resolution for clustering within each sample.
#' @param minibatch Integer specifying the number of columns to process in each minibatch in the SNF analysis.
#' Default is 5000. This option splits the matrix into smaller chunks (minibatch), thus reducing memory usage.
#' @param ncores An integer specifying the number of cores for parallel processing. Default is 1.
#' @param seed An integer used to seed the random number generator.
#'
#' @examples
#' # See https://digitalcytometry.github.io/spatialecotyper/docs/articles/Integration.html
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @import Seurat
#' @import ComplexHeatmap
#'
IntegrateSpatialEcoTyper <- function(SpatialEcoTyper_list,
                                     data_list,
                                     outdir = "./",
                                     normalization.method = "None",
                                     nmf_ranks = 10,
                                     nrun.per.rank = 30,
                                     min.coph = 0.95,
                                     nfeatures = 3000,
                                     min.features = 10,
                                     Region = NULL,
                                     downsample.by.region = TRUE,
                                     subresolution = 30,
                                     minibatch = 5000,
                                     ncores = 1,
                                     seed = 1){

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if(length(SpatialEcoTyper_list)!=length(data_list)){
    stop(" The SpatialEcoTyper_list and data_list should have the same length.")
  }

  if(is.null(names(SpatialEcoTyper_list)) & is.null(names(data_list))){
    warnings("The samples are named as Sample1 -", length(data_list))
    sample_names <- paste0("Sample", 1:length(data_list))
    names(data_list) <- sample_names
    names(SpatialEcoTyper_list) <- sample_names
  }else{
    if(is.null(names(data_list))) names(data_list) = names(SpatialEcoTyper_list)
    if(is.null(names(SpatialEcoTyper_list))) names(SpatialEcoTyper_list) = names(data_list)
    sample_names <- names(data_list)
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

  metadata_list <- mclapply(SpatialEcoTyper_list, function(x){
    obj = FindClusters(x$obj, resolution = subresolution, verbose = FALSE)
    obj$SE = paste0("InitSE", obj$seurat_clusters)
    metadata = AnnotateCells(x$metadata, obj)
    metadata = metadata[!is.na(metadata$SE), ]
    colnames(metadata)[colnames(metadata)=="SE"] = "InitSE"
    metadata
  }, mc.cores = ncores)
  names(metadata_list) <- names(SpatialEcoTyper_list)
  data_list <- lapply(sample_names, function(s){
    idx1 = colnames(data_list[[s]]) %in% rownames(metadata_list[[s]])
    idx2 = colSums(data_list[[s]]>0) >= min.features
    data_list[[s]][, idx1&idx2]
  })
  names(data_list) <- sample_names

  metadata_list <- lapply(sample_names, function(s){
    scmeta = metadata_list[[s]]
    if(!"CID" %in% colnames(scmeta)){
      scmeta = cbind(CID = rownames(scmeta), scmeta)
    }else{
      scmeta$CID = rownames(scmeta)
    }
    scmeta$Sample = s
    scmeta$InitSE = paste0(scmeta$Sample, "..", scmeta$InitSE)
    scmeta[match(colnames(data_list[[s]]), rownames(scmeta)), ]
  })
  ## Merge meta data from all samples
  commoncols <- table(unlist(lapply(metadata_list, colnames)))
  commoncols <- names(commoncols)[commoncols==max(commoncols)]
  commoncols <- unique(c("CID", "Sample", "InitSE", "CellType", commoncols))
  commoncols <- intersect(colnames(metadata_list[[1]]), commoncols)
  commoncols <- setdiff(commoncols, c("SpotID", "Spot.X", "Spot.Y"))
  metadatas <- do.call(rbind, lapply(metadata_list, function(x){ x[, commoncols] }))
  metadatas <- metadatas %>% filter(!is.na(InitSE)) %>% as.data.frame
  names(metadata_list) <- sample_names

  message(Sys.time(), " Construct cell-type-specific gene expression signatures of spatial clusters")
  avgexpr_list <- mclapply(sample_names, function(ss){
    avgexprs = ComputeFCs(data_list[[ss]], metadata_list[[ss]],
                          cluster = "InitSE", Region = Region,
                          scale = TRUE, ncores = ncores)
    return(avgexprs)
  }, mc.cores = ncores)
  #### Select common genes ####
  genes <- table(unlist(lapply(avgexpr_list, rownames)))
  genes <- names(genes)[genes==length(avgexpr_list)]
  genes <- na.omit(genes)
  avgexpr_list <- lapply(avgexpr_list, function(x){ x[genes, ] })
  avgexpr_list <- do.call(cbind, avgexpr_list)

  ## Cluster level meta data
  clustmetas <- metadatas[, -1] %>% group_by(InitSE) %>%
    summarise(across(where(is.character) | where(is.factor) |
                       where(is.logical), mostFrequent)) %>% as.data.frame
  rownames(clustmetas) <- clustmetas$InitSE
  clustmetas <- clustmetas[colnames(avgexpr_list), ]

  message(Sys.time(), " Start integrating SEs across samples ")
  if(!is.null(Region)){
    integrated <- Integrate(avgexpr_list, Region = clustmetas[, Region],
                            downsample.by.region = downsample.by.region,
                            nfeatures = nfeatures, min.features = min.features,
                            minibatch = minibatch, ncores = ncores, seed = seed)
  }else{
    integrated <- Integrate(avgexpr_list, nfeatures = nfeatures,
                            min.features = min.features,
                            minibatch = minibatch, ncores = ncores, seed = seed)
  }
  # saveRDS(integrated, paste0(outdir, "/MultiSE_integrated_priorNMF.rds"))

  message(Sys.time(), " Identify conserved SEs via NMF clustering")
  nmf_res <- nmfClustering(integrated, ranks = nmf_ranks, ncores = ncores,
                           nrun.per.rank = nrun.per.rank,
                           min.coph = min.coph, seed = seed)
  saveRDS(nmf_res, paste0(outdir, "/MultiSE_NMF_results.rds"))
  if(is.list(nmf_res)){
    ggsave(paste0(outdir, "/MultiSE_NMF_Cophenetic_Dynamics.pdf"),
           nmf_res$p, width = 5, height = 4.5)
    ses <- predict(nmf_res$NMFfits[[paste0("K.", nmf_res$bestK)]])
  }else{
    ses <- predict(nmf_res)
  }
  gg <- scale(t(scale(integrated[names(ses), names(ses)])))
  ords <- order.dendrogram(cluster_between_groups(gg, factor(ses)))
  gg <- gg[ords, ords]
  message(Sys.time(), " Save integrated results into ",
          paste0(outdir, "/MultiSE_integrated_final.rds"))
  saveRDS(gg, paste0(outdir, "/MultiSE_integrated_final.rds"))

  ses <- ses[rownames(gg)]
  ses <- factor(ses, levels = unique(ses))
  ann <- data.frame(Sample = gsub("\\.\\..*", "", rownames(gg)),
                    SE = paste0("SE", as.integer(ses)),
                    row.names = rownames(gg))
  if(!is.null(Region)){
    ann$Region = clustmetas[rownames(gg), Region]
    region_cols = getColors(length(unique(ann$Region)), palette = 4)
    names(region_cols) = unique(ann$Region)
  }
  idx <- match(paste0(metadatas$InitSE), rownames(ann))
  metadatas$SE <- ann$SE[idx]
  message(Sys.time(), " Save updated single cell meta data into ",
          paste0(outdir, "/MultiSE_metadata_final.tsv"))
  saveRDS(metadatas, paste0(outdir, "/MultiSE_metadata_final.rds"))
  write.table(metadatas, paste0(outdir, "/MultiSE_metadata_final.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message(Sys.time(), " Visualize the SEs across samples ")
  SE_cols <- getColors(length(unique(ann$SE)), palette = 1)
  names(SE_cols) <- unique(ann$SE)
  sample_cols <- getColors(length(unique(ann$Sample)), palette = 3)
  names(sample_cols) <- unique(ann$Sample)
  if(!is.null(Region)){
    top_ann_col = list(Sample = sample_cols, SE = SE_cols, Region = region_cols)
  }else{
    top_ann_col = list(Sample = sample_cols, SE = SE_cols)
  }

  p <- HeatmapView(gg, breaks = c(0, 0.6, 1.2),
                   colors = c("#ffffd9", "#edf8b1", "#225ea8"),
                   left_ann_col = list(SE = SE_cols), show_left_legend = FALSE,
                   top_ann_col = top_ann_col,
                   show_row_names = FALSE, show_column_names = FALSE,
                   top_ann = ann, left_ann = ann[, 2, drop = FALSE],
                   legend_height = unit(1.5, "cm"))

  pdf(paste0(outdir, "/MultiSE_integrated_final_hmap.pdf"), width = 7, height = 5)
  draw(p)
  drawRectangleAnnotation(ann$SE, ann$SE)
  dev.off()

  ## Spatial landscape of SEs
  p = SpatialView(metadatas, by = "SE")  +
      facet_wrap(~Sample, scales = "free", ncol = 3) +
      scale_color_manual(values = SE_cols)
  nsample = length(unique(metadatas$Sample))
  width = ifelse(nsample>=3, 9, nsample*3)
  height = ceiling(nsample / 3) * 3
  ggsave(paste0(outdir, "/SpatialView_SEs_by_Sample.pdf"), p,
         width = width, height = height)

  ## Cell type composition
  ct_cols <- getColors(length(unique(metadatas$CellType)), palette = 2)
  names(ct_cols) <- unique(metadatas$CellType)
  gg <- metadatas %>% filter(!is.na(SE)) %>% count(SE, CellType, Sample) %>%
    group_by(Sample, SE) %>% mutate(Frac = n / sum(n)) %>%
    ## cell type fractions within each sample
    group_by(SE, CellType) %>% summarise(Frac = mean(Frac))
    ## average cell type fractions across all samples
  p = ggplot(gg, aes(SE, Frac, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = ct_cols) +
    theme_bw(base_size = 14) + coord_flip() +
    labs(y = "Cell type abundance")
  ggsave(paste0(outdir, "/BarView_SEs_CellType_Frac_Avg.pdf"), p,
         width = 6, height = length(unique(metadatas$SE))/4+2)

  ## Mixing among samples
  gg <- metadatas %>% filter(!is.na(SE)) %>% count(SE, Sample) %>%
    group_by(Sample) %>% mutate(n = n / sum(n)) %>%
    group_by(SE) %>% mutate(Frac = n / sum(n))
  ## average cell type fractions across all samples
  p = ggplot(gg, aes(SE, Frac, fill = Sample)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = sample_cols) +
    theme_bw(base_size = 14) + coord_flip() +
    labs(y = "Cell fractions")
  ggsave(paste0(outdir, "/BarView_SEs_Sample_Frac.pdf"), p,
         width = 6, height = length(unique(metadatas$SE))/4+2)

  if(!is.null(Region)){
    metadatas$Region = metadatas[, Region]
    gg <- metadatas %>% filter(!is.na(SE)) %>% count(SE, Region, Sample) %>%
      group_by(Sample, SE) %>% mutate(Frac = n / sum(n)) %>% ## cell type fractions within each sample
      group_by(SE, Region) %>% summarise(Frac = mean(Frac)) ## average cell type fractions across all samples
    p = ggplot(gg, aes(SE, Frac, fill = Region)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = region_cols) +
      theme_bw(base_size = 14) + coord_flip() +
      labs(y = "Fraction")
    ggsave(paste0(outdir, "/BarView_SEs_Region_Frac_Avg.pdf"), p,
           width = 6, height = length(unique(metadatas$SE))/4+2)
  }
}
