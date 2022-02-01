#' @include utilities.R read_alignments.R alignments2long.R
NULL

#' Generate consensus from alignments
#'
#' @description This function takes either Filepaths, named vector alignments, or a data.frame from alignments2long. The output is the consensus sequences in either 'long' or 'stringvector' format.
#' @param input file, character vector or alignments2long output
#' @param outformat 'stringvector' or 'long'
#' @keywords DNA, AA, alignment, consensus
#' @export
#' @examples
#' generate_consensus(input = c(a = 'AAACAAA', a = 'AAACAAA', a = 'AAAAAAA'))

generate_consensus = function(input, outformat = 'stringvector'){
      stopifnot(outformat %in% c('long','stringvector'))

      if(all(file.exists(input))){
	    input = read_alignments(input, naming = 'filenames')
      }
      if(is.character(input)){ 
	    input = alignments2long(input)
      }
      stopifnot(is.list(input) & all(c('Sample','Pos','AA') %in% colnames(input)))

      consens = input %>% dplyr::group_by(Sample,Pos, AA) %>% dplyr::count() %>% dplyr::group_by(Sample,Pos) %>% dplyr::filter(n == max(n)) %>% dplyr::select(-n) %>% dplyr::mutate(nPos = dplyr::n()) 

      #testing for Amino Acid Frequency ties and choosing randomly
      testforties = consens %>% dplyr::group_by(Sample, Pos, nPos) %>% dplyr::summarise() %>% dplyr::ungroup() %>% dplyr::filter(nPos > 1) %>% as.data.frame()
      if(nrow(testforties)>0) {
	    for(x in 1:nrow(testforties)) {
		  warning(paste0('Frequency tie for sample ',testforties[x,1],' position ',testforties[x,2], '. Choosing randomly from the options.'))
	    }
	    consens = consens %>% dplyr::filter(sample(x = c(T,rep(F, n()-1)), n()))
      }
      if(outformat == 'long') return(consens %>% dplyr::select(-nPos))
      if(outformat == 'stringvector') {
	    consens = consens %>% dplyr::select(-nPos) %>% dplyr::group_by(Sample) %>% dplyr::arrange(Pos) %>% dplyr::summarise(Sequence = paste(AA,collapse = ''))
	    outvec = consens$Sequence
	    names(outvec) = consens$Sample
	    return(outvec)
      }
}
