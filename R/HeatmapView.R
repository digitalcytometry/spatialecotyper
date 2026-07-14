#' Draw Heatmap
#'
#' @docType methods
#' @name HeatmapView
#' @rdname HeatmapView
#'
#' @param mat Matrix like object.
#' @param breaks A vector indicating numeric breaks.
#' @param colors A vector of colors which correspond to values in breaks.
#' @param na_col Color for NA values.
#' @param name Heatmap name.
#' @param cluster_rows Same as that in ComplexHeatmap::Heatmap.
#' @param row_dend_side Same as that in ComplexHeatmap::Heatmap.
#' @param cluster_cols Same as that in ComplexHeatmap::Heatmap.
#' @param column_dend_side Same as that in ComplexHeatmap::Heatmap.
#' @param show_row_names Same as that in ComplexHeatmap::Heatmap.
#' @param row_names_side Same as that in ComplexHeatmap::Heatmap.
#' @param row_names_gp Same as that in ComplexHeatmap::Heatmap.
#' @param row_names_rot Same as that in ComplexHeatmap::Heatmap.
#' @param show_column_names Same as that in ComplexHeatmap::Heatmap.
#' @param column_names_side Same as that in ComplexHeatmap::Heatmap.
#' @param column_names_gp Same as that in ComplexHeatmap::Heatmap.
#' @param column_names_rot Same as that in ComplexHeatmap::Heatmap.
#' @param show_legend Whether show annotation legends.
#' @param top_ann A data frame. Each column will be treated as a simple annotation.
#' The data frame must have column names. Can also be a HeatmapAnnotation-class object.
#' @param top_ann_col A list of colors which contain color mapping to df.
#' @param show_top_legend Whether show annotation legends.
#' @param bott_ann Same as top_ann.
#' @param bott_ann_col A list of colors which contain color mapping to df.
#' @param show_bott_legend Whether show annotation legends.
#' @param left_ann Same as top_ann.
#' @param left_ann_col A list of colors which contain color mapping to df.
#' @param show_left_legend Whether show annotation legends.
#' @param right_ann Same as top_ann.
#' @param right_ann_col A list of colors which contain color mapping to df.
#' @param show_right_legend Whether show annotation legends.
#' @param show_ann_name Whether show annotation names.
#' @param annotation_legend_param A list which contains parameters for annotation legends.
#' @param row_split A vector or a data frame by which the rows are split.
#' But if cluster_rows is a clustering object, split can be a single number
#' indicating to split the dendrogram by cutree.
#' @param column_split Same as row_split.
#' @param show_heatmap_legend Whether show legends.
#' @param legend_title Character specifyin the legend title.
#' @param legend_title_position Position of title relative to the legend.
#' topleft, topcenter, leftcenter-rot and lefttop-rot are only for vertical legend and
#' leftcenter, lefttop are only for horizontal legend.
#' @param legend_direction Vertical or horizontal?
#' @param legend_title_gp Graphic parameters of the title.
#' @param legend_labels_gp Graphic parameters for labels.
#' @param legend_height Height of the whole legend body. It is only used for vertical continous legend.
#' @param legend_width Width of the whole legend body. It is only used for horizontal continous legend.
#' @param legend_side Side to put heatmap legend
#' @param ... Other parameters in draw.
#'
#' @return A `Heatmap-class` object.
#' @author Wubing Zhang
#'
#' @examples
#' library(grid)
#' library(SpatialEcoTyper)
#' library(ComplexHeatmap)
#'
#' dat = matrix(rnorm(100), 10)
#' rownames(dat) = letters[1:10]
#' colnames(dat) = letters[11:20]
#' rowann = data.frame(Group = rep(letters[1:2], each=5), index = 1:10)
#' colann = data.frame(Group = rep(letters[1:2], each=5), index = 11:20)
#' HeatmapView(dat, left_ann = rowann, top_ann = colann)
#'
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap draw
#' @importFrom circlize colorRamp2
#' @import grid
#' @export

HeatmapView <- function(mat,
                        breaks = c(0, 0.6, 1.2),
                        colors = c("#ffffd9", "#edf8b1", "#225ea8"),
                        na_col = "grey",
                        name = "hmap",
                        cluster_rows = FALSE,
                        row_dend_side = c("left", "right"),
                        cluster_cols = FALSE,
                        column_dend_side = c("top", "bottom"),
                        show_row_names = TRUE,
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 12),
                        row_names_rot = 0,
                        show_column_names = TRUE,
                        column_names_side = "bottom",
                        column_names_gp = gpar(fontsize = 12),
                        column_names_rot = 90,
                        show_legend = TRUE,
                        top_ann = NULL,
                        top_ann_col = NULL,
                        show_top_legend = show_legend,
                        bott_ann = NULL,
                        bott_ann_col = NULL,
                        show_bott_legend = show_legend,
                        left_ann = NULL,
                        left_ann_col = NULL,
                        show_left_legend = show_legend,
                        right_ann = NULL,
                        right_ann_col = NULL,
                        show_right_legend = show_legend,
                        show_ann_name = TRUE,
                        annotation_legend_param = list(),
                        row_split = NULL,
                        column_split = NULL,
                        show_heatmap_legend = TRUE,
                        legend_title = NULL,
                        legend_title_position = "lefttop",
                        legend_direction = "vertical",
                        legend_title_gp = gpar(fontsize = 12),
                        legend_labels_gp = gpar(fontsize = 12),
                        legend_height = 2,
                        legend_width = 0.3,
                        legend_side = "right",
                        ...){
  mat[which(mat>max(breaks))] = max(breaks)
  mat[which(mat<min(breaks))] = min(breaks)
  colPal = circlize::colorRamp2(breaks, colors)

  if(!(is.null(top_ann)|class(top_ann)=="HeatmapAnnotation")){
    top_ann <- heatmap_annotation(top_ann, palettes = top_ann_col,
                       name = "topann", which = "column",
                       show_legend = show_top_legend,
                       annotation_legend_param = annotation_legend_param,
                       show_annotation_name = show_ann_name)
  }
  if(!(is.null(bott_ann)|class(bott_ann)=="HeatmapAnnotation")){
    bott_ann <- heatmap_annotation(bott_ann, palettes = bott_ann_col,
                       name = "bottann", which = "column",
                       show_legend = show_bott_legend,
                       annotation_legend_param = annotation_legend_param,
                       show_annotation_name = show_ann_name)
  }
  if(!(is.null(left_ann)|class(left_ann)=="HeatmapAnnotation")){
    left_ann <- heatmap_annotation(left_ann, palettes = left_ann_col,
                       name = "leftann", which = "row",
                       show_legend = show_left_legend,
                       annotation_legend_param = annotation_legend_param,
                       show_annotation_name = show_ann_name)
  }
  if(!(is.null(right_ann)|class(right_ann)=="HeatmapAnnotation")){
    right_ann <- heatmap_annotation(right_ann, palettes = right_ann_col,
                       name = "rightann", which = "row",
                       show_legend = show_right_legend,
                       annotation_legend_param = annotation_legend_param,
                       show_annotation_name = show_ann_name)
  }
  tryCatch({
    ht_opt$message = F
  }, error = function(x){})
  ht_opt$verbose = F

  ComplexHeatmap::Heatmap(as.matrix(mat), col = colPal, na_col = na_col, name = name, border = T,
                          cluster_rows = cluster_rows, row_dend_side = row_dend_side,
                          cluster_columns = cluster_cols, column_dend_side = column_dend_side,

                          show_row_names = show_row_names, row_names_side = row_names_side,
                          row_names_gp = row_names_gp, row_names_rot = row_names_rot,
                          show_column_names = show_column_names, column_names_side = column_names_side,
                          column_names_rot = column_names_rot, column_names_gp = column_names_gp,

                          left_annotation = left_ann, right_annotation = right_ann,
                          top_annotation = top_ann, bottom_annotation = bott_ann,

                          row_split = row_split,
                          column_split = column_split,
                          show_heatmap_legend = show_heatmap_legend,
                          heatmap_legend_param = list(
                            title = legend_title,
                            title_gp = legend_title_gp,
                            labels_gp = legend_labels_gp,
                            legend_height = unit(legend_height, "in"),
                            legend_width = unit(legend_width, "in"),
                            title_position = legend_title_position,
                            legend_direction = legend_direction),
                          use_raster = T, raster_device = c("png"), ...)
}

heatmap_annotation <- function(df, palettes = NULL,
                               name = "ann", which = "column",
                               show_legend = TRUE,
                               annotation_legend_param = list(),
                               show_annotation_name = T){
  if(is.null(palettes) | ncol(df)>length(palettes)){
    if(is.null(palettes)) palettes = list()
    columns = setdiff(colnames(df), names(palettes))
    for(col in columns){
      if(all(is.na(df[,col]))) next
      if(is.numeric(df[,col])){
        type_cols <- colorRamp2(quantile(df[,col], c(0.1, 0.5, 0.9)), c("blue", "white", "red"))
      }else{
        type_cols = getColors(length(unique(df[,col])), palette = which(columns==col))
        names(type_cols) = unique(df[,col])
      }
      palettes[[col]] = type_cols
    }
  }
  HeatmapAnnotation(df = df, name = name, col = palettes,
                    show_legend = show_legend, which = which,
                    annotation_legend_param = annotation_legend_param,
                    show_annotation_name = show_annotation_name, border = T)
}

#' Draw Rectangular Annotations for Matching Row/Column Groups in a Heatmap
#'
#' This function overlays rectangular boundaries on a \code{ComplexHeatmap}
#' heatmap to highlight blocks where row and column annotations share the
#' same group label. It detects contiguous segments (blocks) of identical
#' labels in both rows and columns and draws rectangles around all matching
#' row–column block pairs.
#'
#' Unlike simpler implementations, this function correctly handles cases where
#' the same annotation label appears in multiple disjoint segments (e.g., after
#' clustering or reordering).
#' @param ht HeatmapList-class from ComplexHeatmap
#'
#' @param rows A vector of row annotation labels. Length must match the number
#' of rows in the heatmap.
#'
#' @param columns A vector of column annotation labels. Length must match the
#' number of columns in the heatmap.
#'
#' @param col Character. Color of rectangle borders. Default is \code{"black"}.
#'
#' @param heatmap_name Character. Name of the heatmap used in
#' \code{\link[ComplexHeatmap]{decorate_heatmap_body}}. Must match the \code{name}
#' argument used when creating the heatmap. Default is \code{"hmap"}.
#'
#' @param include_na Logical. Whether to treat \code{NA} values as a valid group.
#' If \code{FALSE} (default), \code{NA} values are ignored when drawing rectangles.
#' If \code{TRUE}, \code{NA} is treated as a separate group and rectangles may be
#' drawn for missing labels.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Converts annotation vectors to character format
#'   \item Identifies contiguous blocks of identical labels using run-length encoding
#'   \item Maps block coordinates to the heatmap's normalized coordinate system
#'   \item Draws rectangles around all row–column block pairs sharing the same label
#' }
#'
#' All drawing is performed within \code{\link[ComplexHeatmap]{decorate_heatmap_body}},
#' so the heatmap must already be drawn before calling this function.
#'
#' @return
#' This function is called for its side effects and returns \code{NULL}. It adds
#' graphical elements (rectangles) to an existing heatmap.
#'
#' @examples
#' \dontrun{
#' library(ComplexHeatmap)
#'
#' mat <- matrix(rnorm(100), 10)
#' row_ann <- rep(c("A", "B"), each = 5)
#' col_ann <- rep(c("A", "B"), times = 5)
#'
#' p1 = Heatmap(mat, name = "hmap")
#' p1 = draw(p1)
#' drawRectangleAnnotation(p1,
#'   rows = row_ann,
#'   columns = col_ann,
#'   col = "black"
#' )
#' }
#'
#' @seealso \link[ComplexHeatmap]{decorate_heatmap_body}
#'
#' @importFrom ComplexHeatmap decorate_heatmap_body
#' @importFrom grid grid.rect unit gpar
#'
#' @export
#'
#'
drawRectangleAnnotation <- function(ht, rows, columns, col = "black",
                                    heatmap_name = "hmap",
                                    include_na = FALSE) {

  stopifnot(length(rows) > 0, length(columns) > 0)

  rows <- as.character(rows)
  columns <- as.character(columns)

  if (include_na) {
    rows[is.na(rows)] <- "NA"
    columns[is.na(columns)] <- "NA"
  }

  # ---- helper: contiguous blocks ----
  get_blocks <- function(vec) {
    if (length(vec) == 0) return(NULL)
    r <- rle(vec)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    data.frame(label = r$values, start = starts, end = ends,
               stringsAsFactors = FALSE)
  }

  # ---- label matching ----
  same_label <- function(a, b) {
    if (include_na) {
      identical(a, b)
    } else {
      !is.na(a) && !is.na(b) && a == b
    }
  }

  # Get row and column orders with split info
  ro_list <- ComplexHeatmap::row_order(ht)
  co_list <- ComplexHeatmap::column_order(ht)

  # Ensure they're lists (in case of no splitting)
  if (!is.list(ro_list)) ro_list <- list(ro_list)
  if (!is.list(co_list)) co_list <- list(co_list)

  # Get the number of row and column slices
  n_row_slices <- length(ro_list)
  n_col_slices <- length(co_list)

  # Function to draw rectangles for a specific slice
  draw_for_slice <- function(slice_r, slice_c) {
    # Get indices for this slice
    r_ind <- ro_list[[slice_r]]
    c_ind <- co_list[[slice_c]]

    if (length(r_ind) == 0 || length(c_ind) == 0) return()

    # Get the annotation values for this slice
    rows_slice <- rows[r_ind]
    cols_slice <- columns[c_ind]

    # Find contiguous blocks
    row_blocks <- get_blocks(rows_slice)
    col_blocks <- get_blocks(cols_slice)

    if (is.null(row_blocks) || is.null(col_blocks)) return()

    n_row <- length(rows_slice)
    n_col <- length(cols_slice)

    # Calculate normalized coordinates within the slice
    if (n_row > 0) {
      row_blocks$ymin <- 1 - row_blocks$end / n_row
      row_blocks$ymax <- 1 - (row_blocks$start - 1) / n_row
    }

    if (n_col > 0) {
      col_blocks$xmin <- (col_blocks$start - 1) / n_col
      col_blocks$xmax <- col_blocks$end / n_col
    }

    # Draw rectangles for matching labels
    for (i in seq_len(nrow(row_blocks))) {
      for (j in seq_len(nrow(col_blocks))) {

        if (same_label(row_blocks$label[i], col_blocks$label[j])) {

          # Calculate rectangle position
          x_center <- (col_blocks$xmin[j] + col_blocks$xmax[j]) / 2
          y_center <- (row_blocks$ymin[i] + row_blocks$ymax[i]) / 2
          width <- col_blocks$xmax[j] - col_blocks$xmin[j]
          height <- row_blocks$ymax[i] - row_blocks$ymin[i]

          # Draw the rectangle
          grid::grid.rect(
            x = grid::unit(x_center, "npc"),
            y = grid::unit(y_center, "npc"),
            width = grid::unit(width, "npc"),
            height = grid::unit(height, "npc"),
            gp = grid::gpar(col = col, fill = NA, lwd = 2, lty = 1)
          )
        }
      }
    }
  }

  # Use decorate_heatmap_body for each combination of splits
  for (slice_r in seq_len(n_row_slices)) {
    for (slice_c in seq_len(n_col_slices)) {
      # Call decorate_heatmap_body for each specific slice
      ComplexHeatmap::decorate_heatmap_body(
        heatmap_name,
        {
          draw_for_slice(slice_r, slice_c)
        },
        row_slice = slice_r,
        column_slice = slice_c
      )
    }
  }
}
