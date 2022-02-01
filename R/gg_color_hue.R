#' Color generator function identical to ggplot coloring
#'
#' @description ggplot color function, from https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @param n integer. number of colors
#' @return character vector of color values
#' @keywords color
#' @export



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
