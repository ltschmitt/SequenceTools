#' @include gg_color_hue.R
NULL

#' Switch the order of the colors
#'
#' This function takes a color function and reorders the colors that it outputs.
#' @param n integer. number of colors to produce
#' @param color_function function to use for generating colors
#' @param ordering 'systematic' or 'random'. the method it switches the colors
#' @return character vector with color values
#' @keywords color, plotting
#' @export
#' @examples
#' color_switcher(6, ordering = 'random')

color_switcher = function(n, color_function = gg_color_hue, ordering = "systematic"){
      if(ordering == "systematic") {
	    if(n%%2 == 1){ 
		 ne = n-1
		 nn = c(do.call(c,map2((ceiling(ne/2)+1):ne, 1:ceiling(ne/2) , c)),n)
	    }
	    else nn = do.call(c,map2((ceiling(n/2)+1):n, 1:ceiling(n/2) , c))
      }
      if(ordering == "random") nn = sample(x = 1:n,size = n, replace = F)
      color_function(n)[nn]
}
