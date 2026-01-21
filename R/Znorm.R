#' Weighted / Unweighted Uni-variance Normalization
#'
#' Z-score normalization scales each feature (gene) across all cells to have a mean of 0 and a standard deviation of 1.
#' This function allows for optional weighted univariance normalization, where cells are grouped by a categorical variable.
#'
#' @param mat A matrix of gene expression data, where rows represent genes and columns represent cells.
#' @param groups A character vector specifying the group labels for each cell.
#' If provided, weighted univariance normalization will be performed.
#'
#' @return A matrix of Z-score normalized gene expression data, with rows representing genes and columns representing cells.
#'
#' @examples
#' library(Seurat)
#' library(data.table)
#' library(SpatialEcoTyper)
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("1CoQmU3u8MoVC8RbLUvTDQmOuJJ703HHB"),
#'               "HumanMelanomaPatient1_subset_counts.tsv.gz", overwrite = TRUE)
#' scdata <- fread("HumanMelanomaPatient1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]
#' scdata <- as.matrix(scdata[, -1])
#' tmpobj <- CreateSeuratObject(scdata) %>%
#'           SCTransform(clip.range = c(-10, 10), verbose = FALSE)
#' seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
#' if(seurat_version<5){
#'   normdata <- GetAssayData(tmpobj, "data")
#' }else{
#'   normdata <- tmpobj[["SCT"]]$data
#' }
#' # Z-score normalization
#' znorm_data <- Znorm(normdata)
#' head(znorm_data[, 1:5])
#'
#' # Weighted Z-score normalization
#' drive_download(as_id("12xcZNhpT-xbhcG8kX1QAdTeM9TKeFAUW"),
#'                      "HumanMelanomaPatient1_subset_scmeta.tsv",
#'                     overwrite = TRUE, verbose = FALSE)
#' scmeta <- fread("HumanMelanomaPatient1_subset_scmeta.tsv",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' wtdznorm_data <- Znorm(normdata, groups = scmeta$Region)
#' head(wtdznorm_data[, 1:5])
#'
#' @export
#'
Znorm <- function(mat, groups = NULL){
  if(!(is.matrix(mat)|is(mat, "sparseMatrix"))){
    mat <- as.matrix(mat)
  }
  groups <- as.character(groups)
  rgs <- table(groups)
  if(length(rgs)<2){
    return(Seurat::ScaleData(mat))
  }else{
    wt <- 1 / rgs[groups]
    wt0 <- matrix(rep(wt, nrow(mat)), nrow = nrow(mat), byrow = TRUE)
    wtd.mean <- rowSums(wt0  *  mat) / sum(wt)
    wtd.var <- rowSums(wt0 * (mat - wtd.mean)^2) / (sum(wt) - 1)
    wtd.var[wtd.var<=0] <- 1 ## Prevent 0 variance.
    return((mat - wtd.mean) / sqrt(wtd.var))
  }
}
