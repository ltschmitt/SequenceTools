#' Plot tsne from DNA/AA alignment
#'
#' This is a wrapper that uses a vector of same length strings to calculate a distance matrix which is used to calculate a dimensionality reduction via tSNE. This is then plotted with ggplot
#' @keywords DNA, reverse, complement
#' @export

alignment_tsne = function(alignment, labeling = NULL, distance_metric = 'hamming', perplexity = 30, num_threads = 12){
      distances = stringdist::stringdistmatrix(alignment, method = distance_metric)
      tsne = Rtsne(distances, perplexity = 30, is_distance = T, num_threads = 12)
      dat = tsne$Y %>% data.frame() %>% setNames(c("X", "Y"))
      if(!is.null(labeling)) dat = dat %>% mutate(labeling = labeling)
      
      p = dat %>% ggplot(aes(x = X, y = Y, color = labeling)) + geom_point(alpha = 0.5) + labs(color = '') + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank())
      
      return(p)
}
