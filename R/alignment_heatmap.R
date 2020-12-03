#' @include read_alignments.R utilities.R aa_properties.R
NULL

#' Heatmap of aligned sequences
#'
#' A heatmap of aligned sequences
#' @param x character vector of filenames or named character vector with sequences and labels as names
#' @param nReads integer. number of reads to read in and 
#' @param AAthreshold numeric. AA fraction that indicates the maximum threshold per position to be included in the plot
#' @keywords AA, heatmap
#' @export
#' @examples
#' aa <- c(a = 'YNGTNQ', b = 'YGNMAK')
#' alignment_heatmap(aa)

alignment_heatmap = function(x, nReads, AAthreshold = 0.9 ){
      if(all(file.exists(x))) {
	    alns = read_alignments(input = x, nReads = nReads)
      } else alns = x
      libnames = names(alns)
      names(alns) = paste0('S',1:length(alns))

      # calculate hamming distance and cluster for sorting
      distances = stringdist::stringdistmatrix(alns, method = 'hamming')
      hc = stats::hclust(distances, method = 'ward.D2')
      dend = stats::as.dendrogram(hc) %>% dendextend::set('branches_lwd', 0.5)

      # split sequences into single letters and name the columns according to Positions
      split_alns = do.call(rbind,strsplit(alns,''))
      colnames(split_alns) = paste0(1:ncol(split_alns))
      # reform to longer format and filter by Frequency threshold
      long_alns = split_alns[ stats::order.dendrogram(dend),] %>% as.data.frame() %>% dplyr::mutate(ID = forcats::fct_inorder(as.factor(rownames(.)))) %>% tidyr::pivot_longer(cols = matches('[0-9]+'), names_to = 'Pos', values_to = 'AA') %>% dplyr::mutate(Pos = forcats::fct_inorder(as.factor(Pos))) %>% dplyr::group_by(Pos) %>% dplyr::mutate(Freq = as.vector(prop.table(table(AA))[AA])) %>% dplyr::filter(!any(Freq > AAthreshold))

      #calculate mutual information to get something similar to correlation/similarity
      wide_filtered_alns = long_alns %>% tidyr::pivot_wider(id_cols = 'ID', names_from = 'Pos', values_from = 'AA') %>% dplyr::select(-ID)
      mi_pos = infotheo::mutinformation(wide_filtered_alns)
      # nullify selfcomparisons - they skew the dissimilarity matrix
      null_mat = matrix(c(rep(c(0,rep(1,dim(mi_pos)[1])), dim(mi_pos)[1]-1),0), nrow = dim(mi_pos)[1], ncol = dim(mi_pos)[2])
      mi_pos = mi_pos * null_mat
      # make dissimilarity matrix by substracting with the max value
      diss_mat = stats::as.dist(max(mi_pos)-mi_pos)
      # cluster and prep dendrogram for sorting
      hc_pos = stats::hclust(diss_mat, method = 'ward.D2')
      dend_pos = stats::as.dendrogram(hc_pos) %>% dendextend::set('branches_lwd', 0.5)

      # plot the heatmap
      hmap = long_alns %>% dplyr::group_by(Pos) %>% dplyr::mutate(PosAA = paste0(levels(forcats::fct_infreq(forcats::fct_drop(AA))), collapse = '')) %>% dplyr::ungroup() %>% dplyr::mutate(PosAA = paste0(as.character(Pos),'\t', stringr::str_trunc(PosAA,4, ellipsis = '...')), PosAA = forcats::fct_inorder(as.factor(PosAA)) , PosAA = forcats::fct_relevel(PosAA,levels(PosAA)[stats::order.dendrogram(dend_pos)])) %>% dplyr::ungroup() %>% dplyr::mutate(AA = forcats::fct_drop(AA), AA = forcats::fct_relevel(AA, levels(AA)[order(aa_properties(levels(AA), 'HydropathyIndex' ))])) %>% ggplot2::ggplot(ggplot2::aes(x = ID, y = PosAA, fill = AA)) + ggplot2::geom_raster() + ggplot2::scale_fill_viridis_d(option = 'C') + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(hjust = 0)) + ggplot2::labs(x = 'Sequence Reads', y = 'Residues + AA sorted by Freq')
      # seperate legend for grid plotting
      hmap_legend = cowplot::get_legend( hmap + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)))
      hmap = hmap + ggplot2::theme(legend.position = 'None', plot.margin = ggplot2::unit(c(2,2,2,2),'points'))

      # Sequence alignment dendrogram in ggplot
      gg_hc_alns = ggplot2::ggplot(dendextend::as.ggdend(dend), labels = F, nodes = F, theme = NULL) + ggplot2::theme_void() + ggplot2::scale_x_continuous(expand = c(0,0.01)) + ggplot2::scale_y_continuous(expand = c(0.01,0)) + ggplot2::theme(plot.margin = ggplot2::unit(c(10,0,0,0),'points'))

      # Residue positions mutual information dendrogram in ggplot
      gg_hc_pos = ggplot2::ggplot(dendextend::as.ggdend(dend_pos), labels = F, nodes = F, theme = NULL) + ggplot2::coord_flip(expand = T, clip = 'off') + ggplot2::theme_void()+ ggplot2::scale_x_continuous(expand = c(0.01,0)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::theme(plot.margin = ggplot2::unit(c(0,10,0,2),'points'))

      # null plot for filler spaces can actually be replaced by NULL
      null = ggplot2::ggplot() + ggplot2::theme_void()

      if(length(unique(libnames)) > 1) {
	    # make bar for library association
	    libbar = data.frame(libnames, ID = 1:length(libnames))[stats::order.dendrogram(dend),] %>% dplyr::mutate(ID = forcats::fct_inorder(as.factor(ID))) %>% ggplot2::ggplot(ggplot2::aes(x = ID, y = 1, fill = libnames)) + ggplot2::geom_raster() + ggplot2::theme_void()

	    libbar_legend = cowplot::get_legend( libbar + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)))
	    libbar = libbar + ggplot2::theme(legend.position = 'None', plot.margin = ggplot2::unit(c(2,2,2,2),'points'))
      	    
	    # assemble subplot with library bar
	    p1 = cowplot::plot_grid(null, null, null, null, hmap, gg_hc_pos, nrow = 3, ncol = 2, align = 'h', greedy = T, rel_widths = c(3,1), rel_heights = c(9,1,30))
	    p2 = cowplot::plot_grid(gg_hc_alns, null, libbar, null, hmap, null, nrow = 3, ncol = 2, align = 'v', greedy = T, rel_widths = c(3,1), rel_heights = c(9,1,30))
	    legends = cowplot::plot_grid(libbar_legend, hmap_legend, ncol = 1, rel_heights = c(1,3))
	    cowplot::plot_grid(cowplot::ggdraw(p1) + cowplot::draw_plot(p2), legends, ncol = 2, rel_widths = c(9,1))
      } else {
	    # assemble subplot without library bar
	    p1 = cowplot::plot_grid(null, null, hmap, gg_hc_pos, ncol = 2, align = 'h', greedy = T, rel_widths = c(3,1), rel_heights = c(1,3))
	    p2 = cowplot::plot_grid(gg_hc_alns, null, hmap, null, ncol = 2, align = 'v', greedy = T, rel_widths = c(3,1), rel_heights = c(1,3))
	    # assemble plots and draw on top of each other, hv alignment makes the gap to wide, thats why drawn on top of each other
	    cowplot::plot_grid(cowplot::ggdraw(p1) + cowplot::draw_plot(p2), hmap_legend, ncol = 2, rel_widths = c(9,1))
      }
}
