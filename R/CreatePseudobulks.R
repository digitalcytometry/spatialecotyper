#' Create Pseudo-bulk Mixtures
#'
#' This function generates pseudobulk samples by aggregating single-cell transcriptomics.
#'
#' @param data A matrix of normalized gene expression data (genes x cells). If NULL, counts must be provided.
#' @param groups A named vector indicating the group (e.g., spatial ecotype) for each cell.
#' The names should correspond to the column names of the data matrix.
#' @param counts A matrix of raw counts data (genes x cells). Used to generate normalized data if data is not provided.
#' @param n_mixtures An integer specifying the number of pseudobulk samples to create. Default is 100.
#' @return A list containing two elements:
#' \item{Fracs}{A matrix of the fractions of each group in the pseudobulk samples (rows represent pseudobulk samples, columns represent groups).}
#' \item{Mixtures}{A matrix of pseudobulk gene expression data (genes x pseudobulk samples).}
#' @details
#' \itemize{
#'   \item If `data` is not provided, the function will normalize the `counts` matrix by dividing each column by its sum and multiplying by 10,000.
#'   \item If the maximum value in `data` is less than 80, it assumes the data is in log2 scale and converts it back to non-log scale.
#'   \item The `groups` vector is used to ensure that cells are correctly assigned to their respective groups.
#'   If `groups` does not have names, it is assumed that the names correspond to the column names of the `data` matrix.
#'   \item Pseudobulk samples are created by sampling cells from each group based on predefined fractions,
#'   and then calculating the average expression for each gene in the pseudobulk samples.
#'   \item The pseudobulk data is then normalized using Seurat's `NormalizeData` function.
#' }
#' @examples
#' library(SpatialEcoTyper)
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("15n9zlXed74oeGaO1pythOOM_iWIfuMn2"), "Melanoma_WU2161_counts.rds",
#'                     overwrite = TRUE)
#' counts <- readRDS("Melanoma_WU2161_counts.rds") ## raw counts of scRNA-seq data
#' groups <- sample(paste0("SE", 1:10), ncol(counts), replace = TRUE)
#' names(groups) <- colnames(counts)
#' result = CreatePseudobulks(counts = counts, groups = groups, n_mixtures = 20)
#' head(result$Mixtures[, 1:5]) ## Gene expression matrix of pseudobulks
#' head(result$Fracs) ## SE fractions in pseudobulks
#'
#' @export
#'
CreatePseudobulks <- function(data = NULL, groups, counts = NULL, n_mixtures = 100){
  if(is.null(data)) data = t(t(counts) / colSums(counts) * 10000)
  if(max(data)<80) data = 2^data-1
  if(is.null(names(groups))){
    names(groups) = colnames(data)
  }else{
    groups = groups[colnames(data)]
  }
  fracs <- lapply(sort(unique(groups)), function(x){
    sample(rnorm(n = 10000, mean = 2), n_mixtures)
  })
  fracs <- do.call(cbind, fracs)
  colnames(fracs) <- sort(unique(groups))
  rownames(fracs) <- paste0("Pseudobulk", 1:n_mixtures)
  fracs[fracs<0] <- 0
  fracs <- fracs[rowSums(fracs>0)>2, ]
  fracs <- fracs / rowSums(fracs)
  ncells <- round(fracs*1000)

  pseudobulk <- apply(ncells, 1, function(x){
    ses <- colnames(ncells)
    names(x) <- ses
    cells <- unlist(lapply(ses, function(s){
      sample(names(groups)[groups==s], x[s], replace = TRUE)
    }))
    rowMeans(data[, cells])
  })
  pseudobulk <- Seurat::NormalizeData(pseudobulk, verbose = FALSE)
  return(list(Fracs = fracs, Mixtures = pseudobulk))
}

