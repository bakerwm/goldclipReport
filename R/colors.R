
#' ColorPicker
#'
#' choose colors for barplot
#'
#' @param group integers, pre-defined groups of colors
#' @param n integer, number of colors in each group
#'
#' @return name of colors
#'
#' @import RColorBrewer
#'
#' @export
ColorPicker <- function(group, n = 4) {
  # Pick colors
  #
  # group1: for 3UTR annotation
  # n=9
  d1a <- RColorBrewer::brewer.pal(n = 11, name = "RdGy")
  d1b <- RColorBrewer::brewer.pal(n = 6, name = "RdYlBu")
  g1 <- c(d1a[10:7], d1b[2:6])

  # group2: for Intron
  # n=10
  d2a <- RColorBrewer::brewer.pal(n = 11, name = "PRGn")
  d2b <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
  g2 <- c(d2a[c(1:3, 8:10)], d2b[c(9, 10, 3, 1)])

  # group3: for Intron
  # n=11
  g3 <- c(d2a[c(1:3, 8:9)], d2b[10], d2a[10], d2b[c(9, 1, 3)], "grey40")
  g4 <- c(d2a[c(1:3, 8:9)], d2b[c(1, 3, 10)], d2a[10], d2b[9], "grey40")

  # group5:
  # n=4
  d3a <- RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")
  g5 <- c(d3a[c(4, 3, 1, 2)])
  g6 <- c(d3a[c(2, 3, 5)])

  # group7: paired
  g7 <- RColorBrewer::brewer.pal(n = n, name = "Paired")

  # return values
  va <- list(g1, g2, g3, g4, g5, g6, g7)
  names(va) <- c("anno1", "anno2", "anno3", "anno4", "map", "qc", "others")

  #c <- ifelse(is.null(va[[group]]), g7, c(va[[group]]))
  if(is.null(va[[group]])) {
    return(g7)
  } else {
    return (va[[group]])
  }
}












