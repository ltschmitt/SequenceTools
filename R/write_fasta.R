#' Write sequences as fasta file
#'
#' This function takes a vector of strings, the names for the Headers and the filename to write a fasta file
#' @keywords fasta, sequences, write
#' @export
#' @examples
#' write_fasta(Sequences = c('ACAAACCAGG','AAACCAGA'), 
#' 	Headers = c('Seq1','Seq2'), File = 'sequences.fasta')

write_fasta = function(Sequences, Headers, File){
      Headers = paste0('>',Headers)
      write(paste(Headers,Sequences,sep = "\n"), File)
}
