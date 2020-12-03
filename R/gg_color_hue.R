#' Color generator function identical to ggplot coloring
#'
#' Provide the number of colors you want. The output will be a string of color values
#' @keywords color
#' @export


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
