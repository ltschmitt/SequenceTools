#' Read in sequence alignments
#'
#' This function reads in alignments from fasta files. The sequences need to be aligned, therefore need all have the same length. The function is capable of reading several files, but the contents need to match in string length. The output is a named vector.
#' @param input filename of fasta without linebreaks in sequence
#' @param nReads integer value of number of reads to read in
#' @param naming 'headers' or 'filenames'. Defines where the names come from.
#' @return named character vector
#' @keywords DNA, AA, alignment, read
#' @export


read_alignments = function(input, nReads = -1L, naming = 'headers') {
      stopifnot(all(file.exists(input)))
      stopifnot(naming %in% c('headers','filenames'))

      if(all(dir.exists(input))) input = list.files(input, pattern = '.fasta$|.fa$', full.names = T) else errorCondition('files and directories mixed - not supported')

      message(paste0('reading ',length(input),' sequence file(s)...'))

      sequences = unlist(lapply(1:length(input), function(x) {
	    dat = readLines(input[x], n = nReads*2)
	    seqs = dat[c(F,T)] 
	    if(naming == 'headers') return(stats::setNames(seqs, gsub('>', '', dat[c(T,F)], fixed = T)))
	    rm(dat)
	    if(naming == 'filenames') return(stats::setNames(seqs, rep(gsub('.nt.fasta$|.aa.fasta$|.fasta$|.fa$','',basename(input[x])), length(seqs))))
      }))
      return(sequences)
}
