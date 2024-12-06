% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{getColors}
\alias{getColors}
\title{Generate a List of Colors}
\usage{
getColors(n, palette = 1, categoric = TRUE, exclude = NULL)
}
\arguments{
\item{n}{An integer specifying the number of colors required.}

\item{palette}{An integer specifying the color palette to use. For categorical palettes, valid values are 1 to 7, and for continuous palettes, valid values are 1 to 14.}

\item{categoric}{A logical value indicating whether to use categorical palettes (TRUE) or continuous palettes (FALSE). Default is TRUE.}

\item{exclude}{A character vector of colors to exclude from the generated list. Default is NULL.}
}
\value{
A vector of colors.
}
\description{
This function generates a list of colors based on the specified number and palette type. It supports both categorical and continuous color palettes.
}
\details{
The function uses different sets of color palettes for categorical and continuous data. For categorical data, the function includes palettes such as kelly, cols25, polychrome, glasbey, alphabet2, and alphabet. For continuous data, the function includes palettes like viridis, parula, magma, coolwarm, warmcool, inferno, plasma, and several kovesi linear palettes.

If the number of requested colors exceeds the available colors in the chosen palette, additional colors are sampled from a combined set of all available colors, ensuring the uniqueness of the generated colors.

The function uses `setdiff` to exclude specified colors and ensures that no colors are repeated by setting a seed for reproducibility.
}
\examples{
# Generate 5 categorical colors using the first palette
getColors(5, palette = 1, categoric = TRUE)

# Generate 10 continuous colors using the viridis palette
getColors(10, palette = 1, categoric = FALSE)

}
