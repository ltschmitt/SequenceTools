#' @include utilities.R 
NULL

#' Lookup properties of amino acids
#'
#' @description This function generates a ggplot2 plot that shows an alignment of the provided aligned sequences.
#' @param Alignment named character vector of same length sequences
#' @param Reference named character vector of one sequence
#' @param RemoveConserved logical. removes positions that are all the same in the alignment
#' @param SequenceWrap integer value of positions to show per line
#' @param Coloring 'AA', 'DNA', 'SideChainClass', 'SideChainPolarity', 'SideChainCharge', 'HydropathyIndex', 'MolecularWeight', or 'Occurance'. Defines how the letters are colored.
#' @param SeqRange integers, start and stop position of sequences to use for alignments.
#' @param fontscale numeric value indicating the size of the letters
#' @param order_by_distance logical, sort sequences by hamming distance
#' @return ggplot
#' @keywords AA, amino acid, properties
#' @export
#' @examples
#' aln <- c('ACGTACGT','AAAATTTT','CCCCGGGG')
#' names(aln) <- c('seq1','seq2','seq3')
#' plot_alignment(aln)

plot_alignment = function(Alignment, Reference = NULL, RemoveConserved = T, SequenceWrap = NULL, Coloring = 'AA', SeqRange = NULL, fontscale = 1, order_by_distance = T ) {

      #Colorscheme
      AAcolorscheme = c(
	    'G' = 'darkorange', 'A' = 'darkorange', 'S' = 'darkorange', 'T' = 'darkorange',
	    'C' = 'chartreuse3', 'V' = 'chartreuse3', 'I' = 'chartreuse3', 'L' = 'chartreuse3',
	    'P' = 'chartreuse3', 'F' = 'chartreuse3', 'Y' = 'chartreuse3', 'M' = 'chartreuse3', 
	    'W' = 'chartreuse3',
	    'N' = 'magenta3', 'Q' = 'magenta3', 'H' = 'magenta3',
	    'D' = 'firebrick', 'E' = 'firebrick',
	    'K' = 'dodgerblue3', 'R' = 'dodgerblue3',
	    'B' = 'grey35', 'Z' = 'grey35', 'X' = 'grey35', '-' = 'grey35', '*' = 'grey35')
      DNAcolorscheme = c(
	    'A' = 'darkorange',
	    'C' = 'chartreuse3',
	    'G' = 'firebrick',
	    'T' = 'dodgerblue3',
	    'N' = 'grey35', '-' = 'grey35', '*' = 'grey35')

#      if(all(file.exists(Alignment))) Alignment = ReadAlnFiles(Alignment)
      if(is.character(Alignment)) stopifnot(length(unique(names(Alignment))) == length(Alignment))
      if(is.character(Alignment) & length(Alignment > 1)) Alignment = alignments2long(Alignment, SeqRange)
      stopifnot(is.list(Alignment) & all(c('Sample','Pos','AA') %in% colnames(Alignment)))
      if(length(unique(Alignment$Sample)) > 100) warning('Plotting of over 100 Sequences is not advised') 
      stopifnot(length(unique(Alignment$Sample)) < 200) 
      
      Refname = names(Reference)
      # use reference of alignment if Sample name fits with input
      if(!is.null(Reference)) if(Reference %in% Alignment$Sample) {
	    Refname = Reference
	    Reference = Alignment %>% dplyr::ungroup() %>% dplyr::filter(Sample == Reference) %>% dplyr::arrange(Pos) %>% dplyr::summarise(Sequence = paste(AA,collapse = '')) %>% .$Sequence
      } 
      # use first sequence of alignment if no reference is given
      if(is.null(Reference)) {
	    Refname = unique(Alignment$Sample)[1]
	    Reference = Alignment %>% dplyr::ungroup() %>% dplyr::filter(Sample == unique(Sample)[1]) %>% dplyr::arrange(Pos) %>% dplyr::summarise(Sequence = paste(AA,collapse = '')) %>% .$Sequence
      }

      # assert that reference is one sequence and reference is at least as long as the aligned sequences
      stopifnot(length(Reference) == 1)
      Reference = unlist(strsplit(Reference,'', fixed = T))
      stopifnot(length(Reference) >= max(as.integer((Alignment$Pos))))

      # prepare df for plotting: match to reference according to Pos and replace with dot
      Alignment = Alignment %>% dplyr::ungroup() %>% dplyr::mutate(isref = Sample == Refname, targetAA = Reference[Pos], targetfit = AA == targetAA) 
      if(RemoveConserved & !all(Alignment$targetfit)) Alignment = Alignment %>% dplyr::group_by(Pos) %>% dplyr::filter(!all(targetfit)) %>% dplyr::ungroup()

      Alignment = Alignment %>% dplyr::mutate(label = ifelse(isref, as.character(AA), ifelse(targetfit,'*',as.character(AA)))) %>% dplyr::group_by(Sample) %>% dplyr::mutate(distance = sum(!targetfit)) %>% dplyr::ungroup()
      if(order_by_distance) {Alignment = Alignment %>% dplyr::arrange(isref,-distance) %>% dplyr::mutate(Sample = forcats::fct_inorder(Sample))
      } else Alignment = Alignment %>% dplyr::arrange(isref,dplyr::desc(Sample)) %>% dplyr::mutate(Sample = forcats::fct_inorder(Sample))
      
      #Color
      if(Coloring %in% c('SideChainClass', 'SideChainPolarity', 'SideChainCharge', 'HydropathyIndex', 'MolecularWeight', 'Occurance')){
	    Alignment = Alignment %>% dplyr::mutate(color = AAPropLookup(AA, Coloring))
      } else {
	    Alignment = Alignment %>% dplyr::mutate(color = label)
      }

      #SequenceWrap
      if(!is.null(SequenceWrap)) Alignment = Alignment %>% dplyr::ungroup() %>% dplyr::mutate(Group = findInterval(Pos,Pos[seq(SequenceWrap, SequenceWrap*floor(length(unique(Pos))/SequenceWrap), by = SequenceWrap)]) )

      #xlabels
      xlabels = unique(Alignment$Pos)
      if(!RemoveConserved & length(xlabels)>100) xlabels[xlabels %% 5 != 0 & xlabels != 1] = ''

      Aplot = Alignment %>% dplyr::arrange(Pos) %>% dplyr::mutate(Pos = forcats::fct_inorder(as.factor(Pos))) %>% { ggplot2::ggplot(., ggplot2::aes(y = as.integer(Sample), x = as.integer(Pos), label = label, color = color)) + ggplot2::geom_text(size = 4.5*fontscale) +
      ggplot2::theme(panel.backgroun = ggplot2::element_blank(), legend.position = 'right', axis.text.x.bottom = ggplot2::element_text(angle = -90, vjust = 0.5), axis.text.x.top = ggplot2::element_text(face = 'bold', size = ggplot2::rel(1.3)), axis.ticks.y = ggplot2::element_blank(), axis.ticks.x.top = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = ggplot2::rel(fontscale*0.9)), axis.title = ggplot2::element_text(size = ggplot2::rel(fontscale))) + 
      ggplot2::scale_y_continuous(breaks = 1:length(levels(.$Sample)), labels = levels(.$Sample), name = 'Samples', expand = c(0,0.6), sec.axis = ggplot2::sec_axis(trans = ~., name = 'Distance to Reference', breaks = ggplot2::derive(), labels = 
		  { dplyr::group_by(., Sample, distance) %>% dplyr::summarise() %>% .$distance })) + 
      ggplot2::scale_x_continuous(name = 'Position', breaks = 1:length(levels(.$Pos)), labels = xlabels, expand = c(0,0.6), sec.axis = ggplot2::sec_axis(trans = ~., name = paste0('Reference: ', Refname), breaks = ggplot2::derive(), labels = 
		  { dplyr::group_by(.,Pos,targetAA) %>% dplyr::summarise() %>% .$targetAA }))}

      if(Coloring %in% c('SideChainClass', 'SideChainPolarity', 'SideChainCharge')){
	    Aplot = Aplot + ggplot2::scale_color_brewer(type = 'qual', palette = 'Set1') + ggplot2::labs(color = Coloring)
      }
      if(Coloring %in% c('HydropathyIndex', 'MolecularWeight', 'Occurance')){
	    Aplot = Aplot + ggplot2::scale_color_viridis_c(option = 'plasma') + ggplot2::labs(color = Coloring)
      }
      if(Coloring == 'DNA') Aplot = Aplot + ggplot2::guides(color = ggplot2::guide_legend(ncol = 1)) + ggplot2::scale_color_manual(values = DNAcolorscheme) + ggplot2::labs(color = Coloring)
      if(Coloring == 'AA') Aplot = Aplot + ggplot2::guides(color = ggplot2::guide_legend(ncol = 1)) + ggplot2::scale_color_manual(values = AAcolorscheme) + ggplot2::labs(color = Coloring)

      if(!is.null(SequenceWrap)) Aplot = Aplot + ggplot2::facet_wrap(~Group, scales = 'free_x', ncol = 1) + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank(), panel.spacing.y = ggplot2::unit(2, 'lines'))

      if(any(Alignment$isref)) Aplot = Aplot + ggplot2::theme(axis.text.x.top = ggplot2::element_blank())

      return(Aplot)
}

