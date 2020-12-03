#' Plot tsne from DNA/AA alignment
#'
#' This is a wrapper that uses a vector of same length strings to calculate a distance matrix which is used to calculate a dimensionality reduction via tSNE. This is then plotted with ggplot
#' @param alignment character vector with same length sequences to aligne
#' @param labeling character vector with labels, same order as alignment
#' @param distance_metric see stringdist::stringdist method parameter
#' @param perplexity see Rtsne::Rtsne perplexity parameter
#' @param num_threads see Rtsne::Rtsne num_threads parameter
#' @keywords DNA, reverse, complement
#' @export

alignment_tsne = function(alignment, labeling = NULL, distance_metric = 'hamming', perplexity = 30, num_threads = 12){
      distances = stringdist::stringdistmatrix(alignment, method = distance_metric)
      tsne = Rtsne::Rtsne(distances, perplexity = 30, is_distance = T, num_threads = 12)
      dat = tsne$Y %>% data.frame() %>% stats::setNames(c("X", "Y"))
      if(!is.null(labeling)) dat = dat %>% dplyr::mutate(labeling = labeling)
      
      p = dat %>% ggplot2::ggplot(aes(x = X, y = Y, color = labeling)) + ggplot2::geom_point(alpha = 0.5) + ggplot2::labs(color = '') + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank())
      
      return(p)
}
