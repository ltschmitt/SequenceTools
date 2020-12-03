#' Generate consensus from alignments
#'
#' This function takes either Filepaths, named vector alignments, or a data.frame from alignments2long. The output is the consensus sequences in either 'long' or 'stringvector' format.
#' @keywords DNA, AA, alignment, consensus
#' @export
#' @examples
#' generate_consensus(input = c('alnA.fa','alnB.fa'), nReads = 100, outformat = 'stringvector')

generate_consensus = function(input, outformat = 'stringvector', nReads = 1000){
      stopifnot(outformat %in% c('long','stringvector'))

      if(all(file.exists(input))){
	    input = read_alignments(input, nReads, naming = 'filenames')
      }
      if(is.character(input)){ 
	    # here should be something that reduces the sequences of each Sample to 1000
	    input = alignments2long(input)
      }
      stopifnot(is.list(input) & all(c('Sample','Pos','AA') %in% colnames(input)))

      # here should be something that reduces the sequences of each Sample to 1000

      consens = input %>% group_by(Sample,Pos, AA) %>% count() %>% group_by(Sample,Pos) %>% filter(n == max(n)) %>% select(-n) %>% mutate(nPos = n()) 

      #testing for Amino Acid Frequency ties and choosing randomly
      testforties = consens %>% group_by(Sample, Pos, nPos) %>% summarise() %>% ungroup %>% filter(nPos > 1) %>% as.data.frame
      if(nrow(testforties)>0) {
	    for(x in 1:nrow(testforties)) {
		  warning(paste0('Frequency tie for sample ',testforties[x,1],' position ',testforties[x,2], '. Choosing randomly from the options.'))
	    }
	    consens = consens %>% filter(sample(x = c(T,rep(F, n()-1)), n()))
      }
      if(outformat == 'long') return(consens %>% select(-nPos))
      if(outformat == 'stringvector') {
	    consens = consens %>% select(-nPos) %>% group_by(Sample) %>% arrange(Pos) %>% summarise(Sequence = paste(AA,collapse = ''))
	    outvec = consens$Sequence
	    names(outvec) = consens$Sample
	    return(outvec)
      }
}
