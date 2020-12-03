#' Heatmap of aligned sequences
#'
#' A heatmap of aligned sequences
#' @keywords AA, heatmap
#' @export
#' @examples
#' aa <- c('YNGTNQ', 'YGNMAKA')
#' alignment_heatmap(dna)

alignment_heatmap = function(files, nReads, AAthreshold = 0.9 ){
      alns = read_alignments(input = files, nReads = nReads)
      libnames = names(alns)
      names(alns) = paste0('S',1:length(alns))

      # calculate hamming distance and cluster for sorting
      distances = stringdist::stringdistmatrix(alns, method = 'hamming')
      hc = hclust(distances, method = 'ward.D2')
      dend = as.dendrogram(hc) %>% dendextend::set('branches_lwd', 0.5)

      # split sequences into single letters and name the columns according to Positions
      split_alns = do.call(rbind,strsplit(alns,''))
      colnames(split_alns) = paste0(1:ncol(split_alns))
      # reform to longer format and filter by Frequency threshold
      long_alns = split_alns[ dendextend::order.dendrogram(dend),] %>% as.data.frame() %>% mutate(ID = forcats::fct_inorder(as.factor(rownames(.)))) %>% tidyr::pivot_longer(cols = matches('[0-9]+'), names_to = 'Pos', values_to = 'AA') %>% mutate(Pos = forcats::fct_inorder(as.factor(Pos))) %>% group_by(Pos) %>% mutate(Freq = as.vector(prop.table(table(AA))[AA])) %>% filter(!any(Freq > AAthreshold))

      #calculate mutual information to get something similar to correlation/similarity
      wide_filtered_alns = long_alns %>% tidyr::pivot_wider(id_cols = 'ID', names_from = 'Pos', values_from = 'AA') %>% select(-ID)
      mi_pos = infotheo::mutinformation(wide_filtered_alns)
      # nullify selfcomparisons - they skew the dissimilarity matrix
      null_mat = matrix(c(rep(c(0,rep(1,dim(mi_pos)[1])), dim(mi_pos)[1]-1),0), nrow = dim(mi_pos)[1], ncol = dim(mi_pos)[2])
      mi_pos = mi_pos * null_mat
      # make dissimilarity matrix by substracting with the max value
      diss_mat = as.dist(max(mi_pos)-mi_pos)
      # cluster and prep dendrogram for sorting
      hc_pos = hclust(diss_mat, method = 'ward.D2')
      dend_pos = as.dendrogram(hc_pos) %>% dendextend::set('branches_lwd', 0.5)

      # plot the heatmap
      hmap = long_alns %>% group_by(Pos) %>% mutate(PosAA = paste0(levels(forcats::fct_infreq(forcats::fct_drop(AA))), collapse = '')) %>% ungroup() %>% mutate(PosAA = paste0(as.character(Pos),'\t', stringr::str_trunc(PosAA,4, ellipsis = '…')), PosAA = forcats::fct_inorder(as.factor(PosAA)) , PosAA = forcats::fct_relevel(PosAA,levels(PosAA)[dendextend::order.dendrogram(dend_pos)])) %>% ungroup() %>% mutate(AA = forcats::fct_drop(AA), AA = forcats::fct_relevel(AA, levels(AA)[order(AAPropLookup(levels(AA), 'HydropathyIndex' ))])) %>% ggplot(aes(x = ID, y = PosAA, fill = AA)) + geom_raster() + scale_fill_viridis_d(option = 'C') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(hjust = 0)) + labs(x = 'Sequence Reads', y = 'Residues + AA sorted by Freq')
      # seperate legend for grid plotting
      hmap_legend = get_legend( hmap + guides(fill = guide_legend(ncol = 1)))
      hmap = hmap + theme(legend.position = 'None', plot.margin = unit(c(2,2,2,2),'points'))

      # Sequence alignment dendrogram in ggplot
      gg_hc_alns = ggplot(dendextend::as.ggdend(dend), labels = F, nodes = F, theme = NULL) + theme_void() + scale_x_continuous(expand = c(0,0.01)) + scale_y_continuous(expand = c(0.01,0)) + theme(plot.margin = unit(c(10,0,0,0),'points'))

      # Residue positions mutual information dendrogram in ggplot
      gg_hc_pos = ggplot(dendextend::as.ggdend(dend_pos), labels = F, nodes = F, theme = NULL) + coord_flip(expand = T, clip = 'off') + theme_void()+ scale_x_continuous(expand = c(0.01,0)) + scale_y_continuous(expand = c(0,0)) + theme(plot.margin = unit(c(0,10,0,2),'points'))

      # null plot for filler spaces can actually be replaced by NULL
      null = ggplot()+theme_void()

      if(length(unique(libnames)) > 1) {
	    # make bar for library association
	    libbar = data.frame(libnames, ID = 1:length(libnames))[dendextend::order.dendrogram(dend),] %>% mutate(ID = forcats::fct_inorder(as.factor(ID))) %>% ggplot(aes(x = ID, y = 1, fill = libnames)) + geom_raster() + theme_void()

	    libbar_legend = get_legend( libbar + guides(fill = guide_legend(ncol = 1)))
	    libbar = libbar + theme(legend.position = 'None', plot.margin = unit(c(2,2,2,2),'points'))
      	    
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
	    cowplot::plot_grid(ggdraw(p1) + draw_plot(p2), hmap_legend, ncol = 2, rel_widths = c(9,1))
      }
}
