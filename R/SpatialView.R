#' Visualize Spatial Landscape of Cells / Spots
#'
#' @param obj A Seurat object or a data frame (with cell names as row names).
#' @param by A feature name for plotting, e.g. cell type, region, gene expression.
#' @param X A character specifying the spatial coordinate of x-axis.
#' @param Y A character specifying the the spatial coordinate of y-axis.
#' @param jitter A boolean specifying whether add jitters to the cells.
#' @param slot The slot in Seurat object to pull feature from.
#' @param highlight.cells A vector specifying the cells for highlighting.
#' @param control.cells A vector specifying the control cells as background. If not specified,
#' all the non-highlighting cells will be considered as control.cells.
#' @param bg.downsample An integer specifying the aim for downsampling the control.cells.
#' @param pt.shape Point shape for plotting
#' @param pt.size A numeric specifying the point size of non-control cells.
#' @param pt.alpha A numeric specifying the point transparency of non-control cells.
#' @param bg.color Color of control cells.
#' @param bg.size A numeric specifying the point size of control cells.
#' @param bg.alpha A numeric specifying the point transparency of control cells.
#'
#' @import ggplot2
#'
#' @return A ggplot object.
#'
#' @examples
#' library(data.table)
#' library(Seurat)
#' library(SpatialEcoTyper)
#' library(ggplot2)
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("1CgUOQKrWY_TG61o5aw7J9LZzE20D6NuI"),
#'                     "HumanMelanomaPatient1_subset_scmeta.tsv", overwrite = TRUE)
#' scmeta <- fread("HumanMelanomaPatient1_subset_scmeta.tsv",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#'
#' # Visualize the cell type annotations in the tissue
#' SpatialView(scmeta, by = "CellType", X = "X", Y = "Y") +
#'             scale_color_manual(values = pals::kelly()[-1])
#' SpatialView(scmeta, by = "Region", X = "X", Y = "Y") +
#'             scale_color_brewer(type = "qual", palette = "Set1")
#'
#' @export
#'
#'
SpatialView <- function(obj, by,
                        X = "X", Y = "Y",
                        pt.shape = 20,
                        pt.size = 0.5,
                        pt.alpha = 1,
                        jitter = FALSE,
                        slot = "data",
                        coord.fix = FALSE,
                        highlight.cells = NULL,
                        control.cells = NULL,
                        ncol = 3,
                        bg.downsample = 2000,
                        bg.color = "gray80",
                        bg.size = 0.5,
                        bg.alpha = 0.7){
  data <- obj
  if("Seurat" %in% class(obj)){
    data <- FetchData(obj, na.omit(c(X, Y, by, colnames(obj@meta.data))), slot = slot)
  }

  if(jitter){
    interval <- sort(setdiff(unique(diff(sort(data[, X]))), 0))[1]
    data[, X] <- data[, X] + runif(nrow(data), -interval, interval)
    interval <- sort(setdiff(unique(diff(sort(data[, Y]))), 0))[1]
    data[, Y] <- data[, Y] + runif(nrow(data), -interval, interval)
  }
  colnames(data)[colnames(data)==X] = "X"
  colnames(data)[colnames(data)==Y] = "Y"
  colnames(data)[colnames(data)==by] = "group.by"
  p <- ggplot()
  if(length(highlight.cells)>1){
    highlightdata <- data[rownames(data) %in% highlight.cells, ]
    if(length(control.cells)<2) control.cells <- setdiff(rownames(data), highlight.cells)
    bgdata <- data[rownames(data) %in% control.cells, , drop = FALSE]
    if(nrow(bgdata)>bg.downsample) bgdata <- bgdata[sample(1:nrow(bgdata), bg.downsample), , drop = FALSE]
    p <- p + geom_point(data = bgdata, aes(X, Y), color = bg.color,
                        shape = pt.shape, size = bg.size, alpha = bg.alpha) +
      geom_point(data = highlightdata, aes(X, Y, color = group.by),
                 shape = pt.shape, size = pt.size, alpha = pt.alpha)
  }else if(length(control.cells)>1){
    bgdata <- data[rownames(data) %in% control.cells, , drop = FALSE]
    highlightdata <- data[!rownames(data) %in% control.cells, , drop = FALSE]
    if(nrow(bgdata)>bg.downsample) bgdata <- bgdata[sample(1:nrow(bgdata), bg.downsample), , drop = FALSE]
    p <- p + geom_point(data = bgdata, aes(X, Y), color = bg.color,
                        shape = pt.shape, size = bg.size, alpha = bg.alpha) +
      geom_point(data = highlightdata, aes(X, Y, color = group.by),
                 shape = pt.shape, size = pt.size, alpha = pt.alpha)
  }else{
    p <- ggplot(data, aes(X, Y, color = group.by)) +
      geom_point(size = pt.size, shape = pt.shape, alpha = pt.alpha)
  }
  p <- p + theme_void(base_size = 12)
  if(is.numeric(data[, "group.by"])){
    p <- p + scale_color_gradient(low = "#00204C", high = "#FFE945")
  }else{
    if(length(unique(data[, "group.by"]))<=45){
      colors <- c(pals::kelly()[-(1)], pals::cols25())
      names(colors) <- NULL
      colors <- colors[1:length(unique(data[, "group.by"]))]
      names(colors) <- unique(data[, "group.by"])
      p <- p + scale_color_manual(values = colors) +
        guides(colour = guide_legend(override.aes = list(size=5)))
    }
  }
  if(coord.fix) p <- p + coord_fixed()
  p <- p + labs(color = NULL)
  return(p)
}

