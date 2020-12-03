#' Library distances
#'
#' This function takes a vector of strings with the letters A C G T and converts them to the opposing letter, T G C A. Letters - and N are also accepted but will stay the same. Other letters are not accepted.
#' @keywords DNA, reverse, complement
#' @export

between_lib_distances = function(files, nReads){
      alns = lapply(files, function(x){
	    return(ReadAlnFiles(x, nReads = nReads))
      })

      dists = apply(combn(1:length(files),2),2, function(x) {
	    dists = stringdist::stringdistmatrix(alns[[x[1]]],alns[[x[2]]], method = 'hamming')
	    name = paste(gsub('[$].+$|.aa.fasta','',basename(files[x[1]])), '-', gsub('[$].+$|.aa.fasta','',basename(files[x[2]])))
	    return(data.frame(distances = as.vector(dists), type = name))
      })

      p = do.call(rbind,dists) %>% group_by(type,distances) %>% count %>% group_by(type) %>% mutate(label = ifelse(n == max(n), as.character(type), '')) %>% ggplot(aes(x = distances, y = n, color = type, label = label)) + geom_line() + ggrepel::geom_label_repel()
      return(p)
}

within_lib_distances = function(files, nReads){
      dists = lapply(files, function(x){
	    aln = ReadAlnFiles(x, nReads = nReads)
	    dists = stringdist::stringdistmatrix(aln, method = 'hamming')
	    name = gsub('[$].+$|.aa.fasta','',basename(x))
	    return(data.frame(distances = as.vector(dists), type = name))
      })

      p = do.call(rbind,dists) %>% group_by(type,distances) %>% count %>% ggplot(aes(x = distances, y = n, color = type)) + geom_line()
      return(p)
}
