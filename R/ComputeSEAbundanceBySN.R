#' Neighbor-weight construction and aggregation
#'
#' Build a sparse neighbor-weight matrix via k-NN, optionally filtered by radius
#'
#' @param ref_coords data.frame/matrix of coordinates for the *source* units
#' (rows = e.g. cells or spots).
#' @param query_coords data.frame/matrix of coordinates for the *output*
#'   units (rows = e.g. spatial neighborhoods or spots), columns X, Y. Defaults to `ref_coords`.
#' @param k number of nearest neighbors to search per query point.
#' @param radius numeric distance cutoff; neighbors farther than this are
#'   dropped. Use `Inf` (default) to keep all `k` neighbors regardless of
#'   distance.
#' @param include.self a boolean specifying whether to include self from the knn.
#'
#' @return sparse dgCMatrix, dim = nrow(query_coords) x nrow(ref_coords)
#' @importFrom Matrix Matrix drop0
#' @importFrom RANN nn2
buildKNNWeights <- function(ref_coords, query_coords = ref_coords,
                            k = 10, radius = Inf, include.self = TRUE) {
  query_coords <- as.matrix(query_coords)
  ref_coords   <- as.matrix(ref_coords)
  k <- min(k, nrow(ref_coords))

  knn <- RANN::nn2(data = ref_coords, query = query_coords, k = k)

  weights <- Matrix(0, nrow = nrow(query_coords), ncol = nrow(ref_coords), sparse = TRUE)
  rownames(weights) <- rownames(query_coords)
  colnames(weights) <- rownames(ref_coords)

  for (i in seq_len(nrow(query_coords))) {
    nn_idx  <- knn$nn.idx[i, ]
    nn_dist <- knn$nn.dists[i, ]
    keep <- nn_dist < radius

    if (!include.self) {
      self_pos <- which(nn_idx == i & nn_dist == 0)
      keep[self_pos] <- FALSE
    }

    nn_idx <- nn_idx[keep]
    if (length(nn_idx) == 0) next
    weights[i, nn_idx] <- 1
  }

  drop0(weights)
}

#' Aggregate a values matrix through a neighbor-weight matrix
#'
#' @param cell2se matrix/data.frame, rows = source units (must align with
#'   `colnames(weights)`), columns = features (e.g. SE levels).
#' @param weights sparse matrix, typically from `buildKNNWeights()`.
#' @param sum2one if `TRUE`, normalize each output row to sum to 1 after
#'   aggregation.
#' @param min.cells minimum number of source units (nonzero weight entries)
#'   an output unit must have to be kept in the result.
#'
#' @return dense matrix, dim = (kept output units) x ncol(cell2se)
#' @importFrom Matrix rowSums
aggregateByWeights <- function(cell2se, weights,
                               sum2one = TRUE, min.cells = 5) {
  if(!(ncol(weights)==nrow(cell2se) & all(colnames(weights)%in%rownames(cell2se)))){
    stop("The column names of `weights` do not match row names of `cell2se`.")
  }
  weights = weights[, match(rownames(cell2se), colnames(weights))]
  keep <- Matrix::rowSums(weights > 0) >= min.cells
  if (!any(keep)) {
    stop("No SNs have at least `min.cells` cells.",
         "Try increasing the search radius/k or lowering `min.cells`.")
  }
  weights <- weights[keep, , drop = FALSE]
  weights <- weights / Matrix::rowSums(weights)

  aggr <- weights %*% as.matrix(cell2se)

  if (sum2one) aggr <- aggr / rowSums(aggr)

  as.matrix(aggr)
}


#' Compute SE abundances within spatial neighborhoods (SNs)
#'
#' @param df data.frame, one row per single cell, with at least columns
#'   X, Y (spatial coordinates) and SE.
#' @param spot_coords data.frame of SN centers with columns X, Y
#'   (rownames = SN IDs). If NULL, a regular grid of spatial neighborhoods
#'   is generated automatically, spaced `grid.size` apart.
#' @param radius numeric, radius defining the SN (cells within this distance
#' of a SN center are included).
#' @param grid.size numeric, spacing between adjacent SN centers when `spot_coords` is auto-generated.
#' @param X,Y,SE column names in `df` for x-coord, y-coord, and SE label.
#' @param min.cells minimum number of cells required within a SN
#'   for it to be kept in the output.
#'
#' @return A data frame with the first two columns containing spatial coordinates of SNs,
#' followed by columns containing SE abundances.
#' @importFrom dplyr group_by summarize
#' @export
ComputeSEAbundanceBySN <- function(df, spot_coords = NULL,
                                     radius = 50, grid.size = 50,
                                     k = min(200, nrow(df)),
                                     X = "X", Y = "Y", SE = "SE",
                                     min.cells = 5) {

  # one-hot cell x SE indicator matrix
  se_levels <- sort(unique(df[[SE]]))
  cell2se <- Matrix(0, nrow = nrow(df), ncol = length(se_levels), sparse = TRUE)
  rownames(cell2se) <- rownames(df)
  colnames(cell2se) <- se_levels
  idx <- cbind(seq_len(nrow(df)), match(df[[SE]], se_levels))
  cell2se[idx] <- 1

  sc_coords <- df[, c(X, Y)]

  # auto-generate a regular grid of microregions if none supplied
  if (is.null(spot_coords)) {
    tmpcoord <- as.data.frame(sc_coords)
    colnames(tmpcoord) <- c("X", "Y")
    tmpcoord$SpotID <- paste0("X", round(tmpcoord$X / grid.size),
                              "_Y", round(tmpcoord$Y / grid.size))
    spot_coords <- tmpcoord %>%
      dplyr::group_by(SpotID) %>%
      dplyr::summarize(X = median(X), Y = median(Y)) %>%
      as.data.frame()
    rownames(spot_coords) <- spot_coords$SpotID
  }

  # cells within `radius` of each microregion center
  weights <- buildKNNWeights(query_coords = spot_coords[, c("X", "Y")],
                             ref_coords   = sc_coords[, c(X, Y)],
                             k = k, radius = radius, include.self = TRUE)
  seabunds = aggregateByWeights(cell2se = cell2se, weights = weights,
                                sum2one = TRUE, min.cells = min.cells)
  seabunds = as.data.frame(seabunds)
  seabunds = cbind(spot_coords[match(rownames(seabunds), rownames(spot_coords)),
                               c("X", "Y")], seabunds)
  seabunds
}


#' Smooth spot-level SE abundances across k-nearest-neighbor spots
#'
#' @param se_mat spot x SE matrix (rownames = spot IDs), values indicating
#'   SE abundance/level in each Visium spot.
#' @param spot_coords data.frame/matrix of spot coordinates (rownames = spot
#'   IDs, must overlap with `rownames(se_mat)`), columns given by `X`, `Y`.
#' @param k number of nearest-neighbor spots to average over. Default 7,
#'   matching a Visium hex grid (1 self + 6 immediate neighbors) when
#'   `include.self = TRUE`.
#' @param X,Y column names in `spot_coords` for x/y coordinates.
#' @param include.self if TRUE (default), each spot's own value is included
#'   as one of its k neighbors (so k=7 = self + 6 hex neighbors). If FALSE,
#'   only the k nearest *other* spots are averaged.
#' @param min.neighbors minimum number of neighbor spots required to keep a
#'   spot in the output.
#'
#' @return A data frame with the first two columns containing spatial coordinates of Visium spots,
#' followed by columns containing SE abundances.
#' @export
SmoothSEAbundances <- function(se_mat, spot_coords, k = 7,
                               X = "X", Y = "Y",
                               include.self = TRUE,
                               min.neighbors = 3) {

  spot_coords <- as.data.frame(spot_coords)

  common <- intersect(rownames(se_mat), rownames(spot_coords))
  if (length(common) == 0) {
    stop("rownames do not match in `se_mat` and `spot_coords`.")
  }
  if (length(common) < nrow(se_mat)) {
    warning("Some spots in `se_mat` have no matching coordinates in ",
            "`spot_coords`; those spots will be dropped.")
  }
  se_mat      <- se_mat[common, , drop = FALSE]
  spot_coords <- spot_coords[common, , drop = FALSE]

  weights <- buildKNNWeights(query_coords = spot_coords[, c(X, Y)],
                             ref_coords   = spot_coords[, c(X, Y)],
                             k = k, radius = Inf, include.self = include.self)

  seabunds = aggregateByWeights(cell2se = se_mat, weights = weights,
                                sum2one = FALSE, min.cells = min.neighbors)
  seabunds = as.data.frame(seabunds)
  seabunds = cbind(spot_coords[match(rownames(seabunds), rownames(spot_coords)),
                               c(X, Y)], seabunds)
  seabunds
}

