#' Lookup properties of amino acids
#'
#' This function generates a ggplot2 plot that shows an alignment of the provided aligned sequences. The options:
#' 	Alignment: a string vector containing the aligned sequences
#' 	Reference: provide a sequence that is used as a reference. It will be on top of the plot. If a reference is not provided the first sequence in the alignment is used.
#' 	RemoveConserved: remove the positions in the alignment that are all the same, i.e. conserved.
#' 	SequenceWrap: useful to show long alignments, provide number of positions you want to show per "row"
#' 	Coloring: the way the letters will be colored, options are: 'AA', 'DNA', 'SideChainClass', 'SideChainPolarity', 'SideChainCharge', 'HydropathyIndex', 'MolecularWeight', 'Occurance'
#' 	SeqRange: provide the positions that you want to be plotted
#' 	fontscale: change the size of the letters
#' 	order_by_distance: if true the sequences will be ordered by hamming distance to the reference, if false will provide the sequences according to the order in the alignment vector
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
	    'B' = 'grey35', 'Z' = 'grey35', 'X' = 'grey35', '-' = 'grey35', '•' = 'grey35', '*' = 'grey35')
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
	    Reference = Alignment %>% ungroup() %>% filter(Sample == Reference) %>% arrange(Pos) %>% summarise(Sequence = paste(AA,collapse = '')) %>% .$Sequence
      } 
      # use first sequence of alignment if no reference is given
      if(is.null(Reference)) {
	    Refname = unique(Alignment$Sample)[1]
	    Reference = Alignment %>% ungroup() %>% filter(Sample == unique(Sample)[1]) %>% arrange(Pos) %>% summarise(Sequence = paste(AA,collapse = '')) %>% .$Sequence
      }

      # assert that reference is one sequence and reference is at least as long as the aligned sequences
      stopifnot(length(Reference) == 1)
      Reference = unlist(strsplit(Reference,'', fixed = T))
      stopifnot(length(Reference) >= max(as.integer((Alignment$Pos))))

      # prepare df for plotting: match to reference according to Pos and replace with dot
      Alignment = Alignment %>% ungroup() %>% mutate(isref = Sample == Refname, targetAA = Reference[Pos], targetfit = AA == targetAA) 
      if(RemoveConserved & !all(Alignment$targetfit)) Alignment = Alignment %>% group_by(Pos) %>% filter(!all(targetfit)) %>% ungroup()

      Alignment = Alignment %>% mutate(label = ifelse(isref, as.character(AA), ifelse(targetfit,'*',as.character(AA)))) %>% group_by(Sample) %>% mutate(distance = sum(!targetfit)) %>% ungroup()
      if(order_by_distance) {Alignment = Alignment %>% arrange(isref,-distance) %>% mutate(Sample = forcats::fct_inorder(Sample))
      } else Alignment = Alignment %>% arrange(isref,desc(Sample)) %>% mutate(Sample = forcats::fct_inorder(Sample))
      
      #Color
      if(Coloring %in% c('SideChainClass', 'SideChainPolarity', 'SideChainCharge', 'HydropathyIndex', 'MolecularWeight', 'Occurance')){
	    Alignment = Alignment %>% mutate(color = AAPropLookup(AA, Coloring))
      } else {
	    Alignment = Alignment %>% mutate(color = label)
      }

      #SequenceWrap
      if(!is.null(SequenceWrap)) Alignment = Alignment %>% ungroup() %>% mutate(Group = findInterval(Pos,Pos[seq(SequenceWrap, SequenceWrap*floor(length(unique(Pos))/SequenceWrap), by = SequenceWrap)]) )

      #xlabels
      xlabels = unique(Alignment$Pos)
      if(!RemoveConserved & length(xlabels)>100) xlabels[xlabels %% 5 != 0 & xlabels != 1] = ''

      Aplot = Alignment %>% arrange(Pos) %>% mutate(Pos = forcats::fct_inorder(as.factor(Pos))) %>% { ggplot(., aes(y = as.integer(Sample), x = as.integer(Pos), label = label, color = color)) + geom_text(size = 4.5*fontscale) +
      theme(panel.backgroun = element_blank(), legend.position = 'right', axis.text.x.bottom = element_text(angle = -90, vjust = 0.5), axis.text.x.top = element_text(face = 'bold', size = rel(1.3)), axis.ticks.y = element_blank(), axis.ticks.x.top = element_blank(), axis.text = element_text(size = rel(fontscale*0.9)), axis.title = element_text(size = rel(fontscale))) + 
      scale_y_continuous(breaks = 1:length(levels(.$Sample)), labels = levels(.$Sample), name = 'Samples', expand = c(0,0.6), sec.axis = sec_axis(trans = ~., name = 'Distance to Reference', breaks = derive(), labels = 
		  { group_by(., Sample, distance) %>% summarise %>% .$distance })) + 
      scale_x_continuous(name = 'Position', breaks = 1:length(levels(.$Pos)), labels = xlabels, expand = c(0,0.6), sec.axis = sec_axis(trans = ~., name = paste0('Reference: ', Refname), breaks = derive(), labels = 
		  { group_by(.,Pos,targetAA) %>% summarise() %>% .$targetAA }))}

      if(Coloring %in% c('SideChainClass', 'SideChainPolarity', 'SideChainCharge')){
	    Aplot = Aplot + scale_color_brewer(type = 'qual', palette = 'Set1') + labs(color = Coloring)
      }
      if(Coloring %in% c('HydropathyIndex', 'MolecularWeight', 'Occurance')){
	    Aplot = Aplot + scale_color_viridis_c(option = 'plasma') + labs(color = Coloring)
      }
      if(Coloring == 'DNA') Aplot = Aplot + guides(color = guide_legend(ncol = 1)) + scale_color_manual(values = DNAcolorscheme) + labs(color = Coloring)
      if(Coloring == 'AA') Aplot = Aplot + guides(color = guide_legend(ncol = 1)) + scale_color_manual(values = AAcolorscheme) + labs(color = Coloring)

      if(!is.null(SequenceWrap)) Aplot = Aplot + facet_wrap(~Group, scales = 'free_x', ncol = 1) + theme(strip.background = element_blank(), strip.text = element_blank(), panel.spacing.y = unit(2, 'lines'))

      if(any(Alignment$isref)) Aplot = Aplot + theme(axis.text.x.top = element_blank())

      return(Aplot)
}

