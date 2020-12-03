#' Lookup properties of amino acids
#'
#' This function takes the single letter annotation of amino acids and converts it into the property of choice. Properties available: 'SideChainClass', 'SideChainPolarity', 'SideChainCharge', 'HydropathyIndex', 'MolecularWeight', 'Occurance'.
#' @keywords AA, amino acid, properties
#' @export
#' @examples
#' aa_properties(AA = c('A','A','R','Q','Y'), property = 'HydropathyIndex')

aa_properties = function(AA, property){
      LookupTab = data.frame( row.names = c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"), SideChainClass = c("aliphatic","basic","amide","acid","sulfur-containing","acid","amide","aliphatic","basic aromatic","aliphatic","aliphatic","basic","sulfur-containing","aromatic","cyclic","hydroxyl-containing","hydroxyl-containing","aromatic","aromatic","aliphatic"), SideChainPolarity = c("non polar","basic polar","polar","acidic polar","non polar","acidic polar","polar","non polar","basic polar","non polar","non polar","basic polar","non polar","non polar","non polar","polar","polar","non polar","polar","non polar"), SideChainCharge = c("neutral","positive","neutral","negative","neutral","negative","neutral","neutral","positive(10%)/neutral(90%)","neutral","neutral","positive","neutral","neutral","neutral","neutral","neutral","neutral","neutral","neutral"), HydropathyIndex = c(1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2), MolecularWeight = c(89.094,174.203,132.119,133.104,121.154,147.131,146.146,75.067,155.156,131.175,131.175,146.189,149.208,165.192,115.132,105.093,119.119,204.228,181.191,117.148), Occurance = c(8.76,5.78,3.93,5.49,1.38,6.32,3.90,7.03,2.26,5.49,9.68,5.19,2.32,3.87,5.02,7.14,5.53,1.25,2.91,6.73)/100)
      LookupTab[AA,property]
}
