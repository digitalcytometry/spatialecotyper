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
#' scdata <- fread("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=Melanoma1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]
#' scdata <- as.matrix(scdata[, -1])
#' tmpobj <- CreateSeuratObject(scdata) %>% NormalizeData(verbose = FALSE)
#' normdata <- GetAssayData(tmpobj, layer = "data")
#' # Z-score normalization
#' znorm_data <- Znorm(normdata)
#' head(znorm_data[, 1:5])
#'
#' # Weighted Z-score normalization
#' scmeta <- fread("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=Melanoma1_subset_scmeta.tsv",
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
  if(is.null(groups)){
    return(Seurat::ScaleData(mat))
  }else{
    if(sum(is.na(groups))>0) stop("The group labels can't be NAs.")
    if(length(groups) != ncol(mat)) stop("Length of group labels does not match the number of columns in the mat.")
    groups <- as.character(groups)
    rgs <- table(groups)
    wt <- 1 / rgs[groups]
    wt0 <- matrix(rep(wt, nrow(mat)), nrow = nrow(mat), byrow = TRUE)
    wtd.mean <- rowSums(wt0  *  mat, na.rm = TRUE) / sum(wt)
    wtd.var <- rowSums(wt0 * (mat - wtd.mean)^2, na.rm = TRUE) / (sum(wt) - 1)
    wtd.var[wtd.var<=0] <- 1 ## Prevent 0 variance.
    return((mat - wtd.mean) / sqrt(wtd.var))
  }
}
