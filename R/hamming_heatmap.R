#' @include utilities.R reverse_complement.R read_alignments.R
NULL

#' Heatmap of hamming distances
#'
#' This function makes a heatmap of hamming distances
#' @param x character vector of filenames or named character vector with sequences and labels as names
#' @param nReads integer. number of reads to read in and 
#' @keywords hamming distances
#' @export
#' aa <- c(a = 'YNGTNQ', b = 'YGNMAK')
#' hamming_heatmap(aa)


hamming_heatmap = function(x, nReads){
      if(all(file.exists(x))) {
	    alns = read_alignments(input = x, nReads = nReads)
      } else alns = x
      libnames = names(alns)
      distances = stringdist::stringdistmatrix(alns, method = 'hamming')

      mdistances = as.matrix(distances)
      colnames(mdistances) = paste0('x',1:ncol(mdistances))
      hc = stats::hclust(distances, method = 'ward.D2')
      dend = stats::as.dendrogram(hc) %>% dendextend::set('branches_lwd', 0.5)
      gg_hc = ggplot2::ggplot(dendextend::as.ggdend(dend), labels = F, nodes = F, theme = NULL) + ggplot2::theme_void() + ggplot2::scale_x_continuous(expand = c(0,0.01)) + ggplot2::scale_y_continuous(expand = c(0.01,0)) + ggplot2::theme(plot.margin = ggplot2::unit(c(10,0,0,0),'points'))

      rownames(mdistances) = paste0('Y',1:nrow(mdistances))
      colnames(mdistances) = paste0('X',1:ncol(mdistances))
      mdlong = mdistances[stats::order.dendrogram(dend),stats::order.dendrogram(dend)] %>% as.data.frame() %>% tibble::rownames_to_column('Y') %>% tidyr::pivot_longer(cols = dplyr::starts_with('X'), names_to = 'X') %>% dplyr::mutate(X = forcats::fct_inorder(as.factor(X)), Y = forcats::fct_inorder(as.factor(Y)))

      hmap = mdlong %>% ggplot2::ggplot(ggplot2::aes(x = X, y = Y, fill = value)) + ggplot2::geom_raster() + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + ggplot2::scale_fill_viridis_c(option = 'C') + ggplot2::labs(x = 'Sequence Reads', y = 'Sequence Reads', fill = 'Hamming\nDistance')
      hmap_legend = cowplot::get_legend(hmap)
      hmap = hmap + ggplot2::theme(legend.position = 'None', plot.margin = ggplot2::unit(c(2,2,2,2),'points'))

      if(length(unique(libnames)) > 1) {
	    # make bar for library association
	    libbar = data.frame(libnames, ID = 1:length(libnames))[stats::order.dendrogram(dend),] %>% dplyr::mutate(ID = forcats::fct_inorder(as.factor(ID))) %>% ggplot2::ggplot(ggplot2::aes(x = ID, y = 1, fill = libnames)) + ggplot2::geom_raster() + ggplot2::theme_void() + ggplot2::labs(fill = 'Libraries')
	    libbar_legend = cowplot::get_legend( libbar + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)))
	    libbar = libbar + ggplot2::theme(legend.position = 'None', plot.margin = ggplot2::unit(c(4,4,4,4),'points'))
	    # assemble subplot with library bar
	    cowplot::plot_grid( cowplot::plot_grid(gg_hc, libbar, hmap, ncol = 1, align = 'v', rel_heights = c(9,1,30)) , cowplot::plot_grid(libbar_legend, hmap_legend, ncol = 1, rel_heights = c(1,1)), ncol = 2, rel_widths = c(8,1))
      } else {
	    # assemble subplot without library bar
	    cowplot::plot_grid( cowplot::plot_grid(gg_hc, hmap, ncol = 1, align = 'v', rel_heights = c(1,4)) , hmap_legend, ncol = 2, rel_widths = c(9,1))
      }
}

