#' Convert alignment vector to long table
#'
#' This function takes a named vector of alignment sequences and converts it to a long data.frame with one letter per row. Columns indicate sequence ID and position in sequence. With get_frequencies the output can be condensed to frequency occurances of the letters per position.
#' @keywords DNA, AA, alignment
#' @export
#' @examples
#' alns <- read_alignments(c('alnA.fa','alnB.fa'), nReads = 100, naming = 'filenames')
#' alignments2long(alns)

alignments2long = function(input, seqrange = NULL, get_frequencies = F) {
      output = do.call(rbind,strsplit(input,'',fixed = T)) 
      if(!is.null(seqrange)) output = output[,seqrange]
      output = output %>% reshape2::melt() %>% dplyr::rename(Sample = Var1, Pos = Var2, AA = value) %>% dplyr::mutate(Sample = as.factor(Sample))
      if(get_frequencies) output = output %>% dplyr::group_by(Sample,Pos,AA) %>% dplyr::count() %>% dyplr::group_by(Sample,Pos) %>% dplyr::mutate(Freq = n/sum(n)) %>% ungroup()
      return(output)
}
