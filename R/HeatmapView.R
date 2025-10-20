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
  mat[mat>max(breaks)] = max(breaks)
  mat[mat<min(breaks)] = min(breaks)
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

#' Draw Rectangle Annotations
#'
#' This function draw rectangle grids based on given row and column factors.
#'
#' @param rows A vector of row identifiers.
#' @param columns A vector of column identifiers.
#'
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
#' drawRectangleAnnotation(rowann$Group, colann$Group)
#'
#' @export
drawRectangleAnnotation <- function(rows, columns, col = "black"){
  levels = union(unique(rows), levels(columns))
  rows = factor(as.character(rows), levels = levels)
  columns = factor(as.character(columns), levels = levels)

  dup = sapply(levels, function(lvl) which(rows == lvl)[1]-1)
  fract = dup / length(rows)

  dup = sapply(levels, function(lvl) which(columns == lvl)[1]-1)
  fract2 = dup / length(columns)

  for(i in length(levels):1)
  {

    if(i == 1)
    {
      if(is.na(fract[i]))
      {
        fract[i] = 0
      }
      if(is.na(fract2[i]))
      {
        fract2[i] = 0
      }
    }else{
      if(i == length(levels))
      {
        if(is.na(fract[i]))
        {
          fract[i] = 1
        }
        if(is.na(fract2[i]))
        {
          fract2[i] = 1
        }
      }else{
        if(is.na(fract[i]))
        {
          fract[i] = fract[i+1]
        }
        if(is.na(fract2[i]))
        {
          fract2[i] = fract2[i+1]
        }
      }
    }
  }
  height =  c(fract[-1], 1) - fract
  width =  c(fract2[-1], 1) - fract2
  rect <- list(y=c(1-fract, 0), h=height, x=c(fract2,1), w=width)
  ## Draw grids
  for(i in 1:(length(rect$x)-1)){
    decorate_heatmap_body("hmap", {
      grid.lines(x = unit(rep(rect$x[i], 2), "native"), y = unit(rect$y[i:(i+1)], "native"),
                 gp = gpar(col = col, lty = 1, lwd = 2))
      grid.lines(x = unit(rep(rect$x[i+1], 2), "native"), y = unit(rect$y[i:(i+1)], "native"),
                 gp = gpar(col = col, lty = 1, lwd = 2))
      grid.lines(x = unit(rect$x[i:(i+1)], "native"), y = unit(rep(rect$y[i],2), "native"),
                 gp = gpar(col = col, lty = 1, lwd = 2))
      grid.lines(x = unit(rect$x[i:(i+1)], "native"), y = unit(rep(rect$y[i+1], 2), "native"),
                 gp = gpar(col = col, lty = 1, lwd = 2))
    })
  }
}

