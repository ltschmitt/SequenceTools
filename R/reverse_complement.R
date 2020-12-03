#' Make the reverse complement of strings of DNA
#'
#' This function takes a vector of strings with the letters A C G T and converts them to the opposing letter, T G C A. Letters - and N are also accepted but will stay the same. Other letters are not accepted.
#' @keywords DNA, reverse, complement
#' @export
#' @examples
#' dna <- c('ACGT', 'CGCGAAA')
#' reverse_complement(dna)


reverse_complement <- function(x, reverse = T, collapsed = T){
	bases=c("A","C","G","T",'-','N')
	comp_bases=c('T','G','C','A','-','N')
	split_sequence = strsplit(toupper(x),NULL)
	out = lapply(split_sequence,function(i){ 
		processed = comp_bases[match(i, bases)]
		if(reverse == T & collapsed == T) {processed = rev(processed)}
		paste0(processed,collapse = '') })
	if(reverse == T & collapsed == F) out = rev(out)
	output = do.call(c, out)
	return(output)
}
