#' Identify the most frequent category in a vector
#'
#' This function takes a vector as input and returns the most frequent
#' category or value. It works by converting the input into a table of
#' frequencies, sorting the frequencies in descending order, and returning
#' the most common value as a character string.
#'
#' @param x A vector of values (e.g., categorical or character data).
#'
#' @return A character string representing the most frequent value in the input vector.
#' @export
#'
mostFrequent <- function(x){ # compute the most frequent category
  table(as.character(x)) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
    slice(1) %>% pull(1) %>% as.character
}

#' Transform a Sparse Matrix to Rank Space (Rank Non-zeros in Each Column)
#'
#' This function ranks the non-zero elements of a sparse matrix within each column.
#' The ranks are scaled by dividing them by the number of rows in the matrix,
#' resulting in ranks normalized between 0 and 1. The original zero values are left unchanged.
#'
#' @param sparseMat A sparse matrix
#'
#' @return A sparse matrix of the same dimensions as the input, where the non-zero elements
#'         are replaced by their ranks within each column, normalized by the number of rows.
#'
#' @import Matrix
#'
#' @examples
#' library(Matrix)
#' # Create a sample sparse matrix
#' sparseMat <- rsparsematrix(5, 5, density = 0.5)
#' # Apply the ranking function
#' rankedSparseMat <- rankSparse(sparseMat)
#' print(rankedSparseMat)
#'
#' @export
#'
rankSparse <- function(sparseMat){
  if(is.matrix(sparseMat)){
    sparseMat <- apply(sparseMat, 2, rank) / nrow(sparseMat)
  }else{
    non_zero_indices <- which(sparseMat != 0, arr.ind = TRUE)
    non_zero_values <- sparseMat@x

    # Step 2: Rank non-zero elements within each column
    sparseMat@x <- ave(non_zero_values, non_zero_indices[, "col"], FUN = rank)
    sparseMat = sparseMat / nrow(sparseMat)
  }
  return(sparseMat)
}



#' Generate a List of Colors
#'
#' This function generates a list of colors based on the specified number and palette type. It supports both categorical and continuous color palettes.
#'
#' @param n An integer specifying the number of colors required.
#' @param palette An integer specifying the color palette to use. For categorical palettes, valid values are 1 to 7, and for continuous palettes, valid values are 1 to 14.
#' @param categoric A logical value indicating whether to use categorical palettes (TRUE) or continuous palettes (FALSE). Default is TRUE.
#' @param exclude A character vector of colors to exclude from the generated list. Default is NULL.
#'
#' @return A vector of colors.
#'
#' @details
#' The function uses different sets of color palettes for categorical and continuous data. For categorical data, the function includes palettes such as kelly, cols25, polychrome, glasbey, alphabet2, and alphabet. For continuous data, the function includes palettes like viridis, parula, magma, coolwarm, warmcool, inferno, plasma, and several kovesi linear palettes.
#'
#' If the number of requested colors exceeds the available colors in the chosen palette, additional colors are sampled from a combined set of all available colors, ensuring the uniqueness of the generated colors.
#'
#' The function uses `setdiff` to exclude specified colors and ensures that no colors are repeated by setting a seed for reproducibility.
#'
#' @examples
#' # Generate 5 categorical colors using the first palette
#' getColors(5, palette = 1, categoric = TRUE)
#'
#' # Generate 10 continuous colors using the viridis palette
#' getColors(10, palette = 1, categoric = FALSE)
#'
#' @export
#'
getColors <- function(n, palette = 1, categoric = TRUE, exclude = NULL){
  require(pals)
  if(categoric){
    allcolors <- c(kelly()[c(3:6,2,7:22, 1)], cols25(), polychrome(), glasbey(), alphabet2(), alphabet())
    if(palette==1){
      colors <- kelly()[c(3:6,2,7:22, 1)]
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(1)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else if(palette==2){
      colors <- cols25()
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(2)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else if(palette==3){
      colors <- polychrome()
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(3)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else if(palette==4){
      colors <- brewer.set1(8)
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(4)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else if(palette==5){
      colors <- brewer.set2(8)
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(5)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else if(palette==6){
      colors <- brewer.dark2(8)
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(6)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else if(palette==7){
      colors <- brewer.dark2(12)
      colors <- setdiff(colors, exclude)
      if(n>length(colors)){
        set.seed(7)
        colors <- c(colors, sample(setdiff(allcolors, c(colors,exclude)), n-length(colors), replace = TRUE))
      }else colors <- colors[1:n]
    }else{
      colors <- allcolors
    }
    names(colors) <- NULL
    colors[1:n]
  }else{
    if(palette==1){
      viridis(n)
    }else if(palette==2){
      parula(n)
    }else if(palette==3){
      magma(n)
    }else if(palette==4){
      coolwarm(n)
    }else if(palette==5){
      warmcool(n)
    }else if(palette==6){
      inferno(n)
    }else if(palette==7){
      plasma(n)
    }else if(palette==8){
      kovesi.linear_bgy_10_95_c74(n)
    }else if(palette==9){
      kovesi.linear_bgyw_15_100_c67(n)
    }else if(palette==10){
      kovesi.linear_blue_95_50_c20(n)
    }else if(palette==11){
      kovesi.linear_bmw_5_95_c86(n)
    }else if(palette==12){
      kovesi.linear_bmy_10_95_c78(n)
    }else if(palette==13){
      kovesi.linear_kry_5_95_c72(n)
    }else if(palette==14){
      kovesi.linear_bgyw_15_100_c67(n)
    }else{
      stop("Palette not found")
    }
  }
}
