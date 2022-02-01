#' Read in sequence alignments
#'
#' This function reads in alignments from fasta files. The sequences need to be aligned, therefore need all have the same length. The function is capable of reading several files, but the contents need to match in string length. The output is a named vector.
#' @param input filename of fasta without linebreaks in sequence
#' @param nReads integer value of number of reads to read in
#' @param naming 'headers', 'filenames' or 'filenames-header'. Defines where the names come from.
#' @return named character vector
#' @keywords DNA, AA, alignment, read
#' @export

read_alignments = function(input, nReads = -1L, naming = 'headers', fileformat = 'singleline') {
      stopifnot(all(file.exists(input)))
      stopifnot(naming %in% c('headers','filenames','filenames-headers'))
      stopifnot(fileformat %in% c('singleline','multiline'))

      if(all(dir.exists(input))) input = list.files(input, pattern = '.fasta$|.fa$', full.names = T) else errorCondition('files and directories mixed - not supported')

      message(paste0('reading ',length(input),' sequence file(s)...'))

      sequences = unlist(lapply(1:length(input), function(x) {
	    dat = readLines(input[x], n = nReads*2)
	    dat = dat[dat != '']
	    if(fileformat == 'singleline'){
		  seqs = dat[c(F,T)] 
		  if(naming == 'headers') return(stats::setNames(seqs, gsub('>', '', dat[c(T,F)], fixed = T)))
		  if(naming == 'filenames-headers') return(stats::setNames(seqs, paste0(rep(gsub('.nt.fasta$|.aa.fasta$|.fasta$|.fa$','',basename(input[x])), length(seqs)), gsub('>', '', dat[c(T,F)], fixed = T), sep = '--')))
		  rm(dat)
		  if(naming == 'filenames') return(stats::setNames(seqs, rep(gsub('.nt.fasta$|.aa.fasta$|.fasta$|.fa$','',basename(input[x])), length(seqs))))
	    } else {
		  seqs = grep('^[^>]',dat)
		  seqs = tibble::tibble(Group = cumsum(c(1,diff(seqs)) !=1), Sequence = dat[seqs]) %>% dplyr::group_by(Group) %>% dplyr::summarise(Sequence = paste0(Sequence,collapse = '')) %>% dplyr::pull(Sequence)
		  if(naming == 'headers') return(stats::setNames(seqs, gsub('>', '', grep('^>',dat,value=T), fixed = T)))
		  if(naming == 'filenames-headers') return(stats::setNames(seqs, paste0(rep(gsub('.nt.fasta$|.aa.fasta$|.fasta$|.fa$','',basename(input[x])), length(seqs)), gsub('>', '', grep('^>',dat, value = T), fixed = T), sep = '--')))
		  rm(dat)
		  if(naming == 'filenames') return(stats::setNames(seqs, rep(gsub('.nt.fasta$|.aa.fasta$|.fasta$|.fa$','',basename(input[x])), length(seqs))))
	    }
      }))
      return(sequences)
}
