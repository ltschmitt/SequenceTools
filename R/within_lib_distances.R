#' @include utilities.R
NULL

#' Within library distances
#'
#' Plot the all vs all hamming distances of alignments from one or more libraries
#' @param alignments character vector of sequences, must be same length
#' @param labeling character vector of library association
#' @keywords Library, comparison
#' @export

within_lib_distances = function(alignments, labeling) {
      dists = lapply(unique(labeling), function(x) {
	    alns = alignments[labeling == x]
	    dists = stringdist::stringdistmatrix(alns, method = "hamming")
	    return(data.frame(distances = as.vector(dists), type = x))
      })
      p = do.call(rbind, dists) %>% dplyr::group_by(type, distances) %>% dplyr::count() %>% ggplot2::ggplot(ggplot2::aes(x = distances, y = n, color = type)) + ggplot2::geom_line()
      return(p)
}
