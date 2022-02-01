#' @include utilities.R
NULL

#' Between library distances
#'
#' Plot the hamming distances between alignments from two or more libraries
#' @param alignments character vector of sequences, must be same length
#' @param labeling character vector of library association
#' @keywords Library, comparison
#' @export

between_lib_distances = function(alignments, labeling){
      ulabels = unique(labeling)
      alns = lapply(unique(labeling), function(x){
	    return(alignments[labeling == x])
      })

      dists = apply(utils::combn(1:length(ulabels),2),2, function(x) {
	    dists = stringdist::stringdistmatrix(alns[[x[1]]],alns[[x[2]]], method = 'hamming')
	    name = paste(ulabels[x[1]], '-', ulabels[x[2]])
	    return(data.frame(distances = as.vector(dists), type = name))
      })

      p = do.call(rbind,dists) %>% dplyr::group_by(type,distances) %>% dplyr::count() %>% dplyr::group_by(type) %>% dplyr::mutate(label = ifelse(n == max(n), as.character(type), '')) %>% ggplot2::ggplot(ggplot2::aes(x = distances, y = n, color = type, label = label)) + ggplot2::geom_line() + ggrepel::geom_label_repel()
      return(p)
}
