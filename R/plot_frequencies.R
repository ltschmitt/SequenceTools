#' @include utilities.R alignments2long.R
NULL

#' Plot Frequencies of libraries
#'
#' Plot the positional frequencies of amino acids in sequence alignments and color it according to match with reference
#' @param alignment character vector with same length sequences to aligne
#' @param reference character vector with one sequence acting as reference to compare the AA to
#' @keywords DNA, reverse, complement
#' @export

plot_frequencies = function(alignment, reference){
      seqfreq = alignments2long(alignment, get_frequencies = T)
      refsplit = unlist(strsplit(unname(reference),''))

      p = seqfreq %>% mutate(same_as_ref = ifelse(refsplit[Pos] == AA, 'yes', 'no')) %>% group_by(Sample, Pos) %>% ungroup() %>% arrange(Pos,desc(same_as_ref)) %>% mutate(Pos = fct_inorder(as.factor(Pos)), same_as_ref = fct_inorder(as.factor(same_as_ref)), label = ifelse(Freq > 0.1, AA, '')) %>% ggplot(aes(x = Pos, y = Freq, fill = same_as_ref, label = label)) + geom_col( color = 'black') + geom_text(position = 'stack', vjust = 1.2) + scale_fill_manual(values = rev(c('gold', 'steelblue3'))) + theme(axis.text.x = element_text(angle = 90))
      return(p)
}
