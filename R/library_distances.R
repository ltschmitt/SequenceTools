#' @include read_alignments.R utilities.R
NULL

#' Library distances
#'
#' This function takes a vector of strings with the letters A C G T and converts them to the opposing letter, T G C A. Letters - and N are also accepted but will stay the same. Other letters are not accepted.
#' @param files character vector of filenames
#' @param nReads integer. number of reads to read in and 
#' @keywords DNA, reverse, complement
#' @export

between_lib_distances = function(files, nReads){
      alns = lapply(files, function(x){
	    return(read_alignments(x, nReads = nReads))
      })

      dists = apply(utils::combn(1:length(files),2),2, function(x) {
	    dists = stringdist::stringdistmatrix(alns[[x[1]]],alns[[x[2]]], method = 'hamming')
	    name = paste(gsub('[$].+$|.aa.fasta','',basename(files[x[1]])), '-', gsub('[$].+$|.aa.fasta','',basename(files[x[2]])))
	    return(data.frame(distances = as.vector(dists), type = name))
      })

      p = do.call(rbind,dists) %>% dplyr::group_by(type,distances) %>% dplyr::count %>% dplyr::group_by(type) %>% dplyr::mutate(label = ifelse(n == max(n), as.character(type), '')) %>% ggplot2::ggplot(ggplot2::aes(x = distances, y = n, color = type, label = label)) + ggplot2::geom_line() + ggrepel::geom_label_repel()
      return(p)
}

within_lib_distances = function(files, nReads){
      dists = lapply(files, function(x){ aln = read_alignments(x, nReads = nReads)
	    dists = stringdist::stringdistmatrix(aln, method = 'hamming')
	    name = gsub('[$].+$|.aa.fasta','',basename(x))
	    return(data.frame(distances = as.vector(dists), type = name))
      })

      p = do.call(rbind,dists) %>% dplyr::group_by(type,distances) %>% dplyr::count %>% ggplot2::ggplot(ggplot2::aes(x = distances, y = n, color = type)) + ggplot2::geom_line()
      return(p)
}
