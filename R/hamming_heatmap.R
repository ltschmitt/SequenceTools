#' Heatmap of hamming distances
#'
#' This function makes a heatmap of hamming distances
#' @keywords hamming distances
#' @export


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

hamming_heatmap = function(files, nReads){
      alns = ReadAlnFiles(input = files, nReads = nReads)
      libnames = names(alns)
      distances = stringdist::stringdistmatrix(alns, method = 'hamming')

      mdistances = as.matrix(distances)
      colnames(mdistances) = paste0('x',1:ncol(mdistances))
      hc = hclust(distances, method = 'ward.D2')
      dend = as.dendrogram(hc) %>% dendextend::set('branches_lwd', 0.5)
      gg_hc = ggplot(dendextend::as.ggdend(dend), labels = F, nodes = F, theme = NULL) + theme_void() + scale_x_continuous(expand = c(0,0.01)) + scale_y_continuous(expand = c(0.01,0)) + theme(plot.margin = unit(c(10,0,0,0),'points'))

      rownames(mdistances) = paste0('Y',1:nrow(mdistances))
      colnames(mdistances) = paste0('X',1:ncol(mdistances))
      mdlong = mdistances[dendextend::order.dendrogram(dend),dendextend::order.dendrogram(dend)] %>% as.data.frame() %>% rownames_to_column('Y') %>% tidyr::pivot_longer(cols = starts_with('X'), names_to = 'X') %>% mutate(X = forcats::fct_inorder(as.factor(X)), Y = forcats::fct_inorder(as.factor(Y)))

      hmap = mdlong %>% ggplot(aes(x = X, y = Y, fill = value)) + geom_raster() + theme(axis.text = element_blank(), axis.ticks = element_blank()) + scale_fill_viridis_c(option = 'C') + labs(x = 'Sequence Reads', y = 'Sequence Reads', fill = 'Hamming\nDistance')
      hmap_legend = get_legend(hmap)
      hmap = hmap + theme(legend.position = 'None', plot.margin = unit(c(2,2,2,2),'points'))

      if(length(unique(libnames)) > 1) {
	    # make bar for library association
	    libbar = data.frame(libnames, ID = 1:length(libnames))[dendextend::order.dendrogram(dend),] %>% mutate(ID = forcats::fct_inorder(as.factor(ID))) %>% ggplot(aes(x = ID, y = 1, fill = libnames)) + geom_raster() + theme_void() + labs(fill = 'Libraries')
	    libbar_legend = get_legend( libbar + guides(fill = guide_legend(ncol = 1)))
	    libbar = libbar + theme(legend.position = 'None', plot.margin = unit(c(4,4,4,4),'points'))
	    # assemble subplot with library bar
	    cowplot::plot_grid( cowplot::plot_grid(gg_hc, libbar, hmap, ncol = 1, align = 'v', rel_heights = c(9,1,30)) , cowplot::plot_grid(libbar_legend, hmap_legend, ncol = 1, rel_heights = c(1,1)), ncol = 2, rel_widths = c(8,1))
      } else {
	    # assemble subplot without library bar
	    cowplot::plot_grid( cowplot::plot_grid(gg_hc, hmap, ncol = 1, align = 'v', rel_heights = c(1,4)) , hmap_legend, ncol = 2, rel_widths = c(9,1))
      }
}

