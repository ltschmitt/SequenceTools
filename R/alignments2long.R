#' @include utilities.R
NULL

#' Convert alignment vector to long table
#'
#' This function takes a named vector of alignment sequences and converts it to a long data.frame with one letter per row. Columns indicate sequence ID and position in sequence. With get_frequencies the output can be condensed to frequency occurances of the letters per position.
#' @param input named character vector
#' @param seqrange integer vector describing the positions in the sequences to use
#' @param get_frequencies logical. calculat frequencies or not
#' @return tibble with colnames c('Sample','Pos','AA') or c('Sample','Pos','AA','Freq')
#' @keywords DNA, AA, alignment
#' @export
#' @examples
#' alignments2long(c(a = 'AAACCCGT', b = 'AAAGGGCT', c = 'AGAGAGAG'))

alignments2long = function(input, seqrange = NULL, get_frequencies = F) {
      output = do.call(rbind,strsplit(input,'',fixed = T)) 
      if(!is.null(seqrange)) output = output[,seqrange]
      output = output %>% reshape2::melt() %>% dplyr::rename(Sample = Var1, Pos = Var2, AA = value) %>% dplyr::mutate(Sample = as.factor(Sample))
      if(get_frequencies) output = output %>% dplyr::group_by(Sample,Pos,AA) %>% dplyr::count() %>% dplyr::group_by(Sample,Pos) %>% dplyr::mutate(Freq = n/sum(n)) %>% dplyr::ungroup()
      return(output)
}
