#' Partition Tissue into Subregions
#'
#' This function partitions a tissue dataset into smaller subregions based on spatial coordinates.
#' The partitioning is determined by dividing the range of the X and Y coordinates into
#' specified numbers of rows (`nrow`) and columns (`ncol`). Each subregion is assigned a
#' unique identifier, facilitating independent analysis of smaller subsets of the tissue.
#'
#' @param meta A data frame containing metadata for the tissue samples. It must include
#'        columns representing spatial coordinates (`X` and `Y`).
#' @param nrow An integer specifying the number of rows to divide the tissue into. Default is 2.
#' @param ncol An integer specifying the number of columns to divide the tissue into. Default is 2.
#' @param X A string representing the name of the column containing X coordinates. Default is "X".
#' @param Y A string representing the name of the column containing Y coordinates. Default is "Y".
#'
#' @return A modified data frame with an additional column, `Partition`, that labels each row
#'         with the corresponding subregion identifier in the format "row_column".
#'
#' @examples
#' # Example metadata with spatial coordinates
#' meta <- data.frame(
#'   X = c(10, 20, 30, 40, 50),
#'   Y = c(15, 25, 35, 45, 55)
#' )
#'
#' # Partition the tissue into 2 rows and 2 columns
#' PartitionTissue(meta, nrow = 2, ncol = 2)
#'
#' @keywords internal
PartitionTissue <- function(meta, nrow = 2, ncol = 2, X = "X", Y = "Y"){
  xlen = max(meta[, X]) - min(meta[, X])
  ylen = max(meta[, Y]) - min(meta[, Y])
  row_width = xlen / nrow
  col_width = ylen / ncol
  row_group = floor((meta[, X] - min(meta[, X])) / row_width)
  col_group = floor((meta[, Y] - min(meta[, Y])) / col_width)
  meta$Partition = paste0(row_group, "_", col_group)
  meta
}
