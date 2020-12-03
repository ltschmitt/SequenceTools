#' Write sequences as fasta file
#'
#' This function takes a vector of strings, the names for the Headers and the filename to write a fasta file
#' @param Sequences character vector containg the sequences to write into the file
#' @param Headers character vector containg the header names to write into the file, same order as sequences
#' @param File Filename to write to
#' @keywords fasta, sequences, write
#' @export

write_fasta = function(Sequences, Headers, File){
      Headers = paste0('>',Headers)
      write(paste(Headers,Sequences,sep = "\n"), File)
}
