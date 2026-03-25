##' Get microRNA-target Interactions
##'
##' Retrieve predicted and validated miRNA-target interactions using TargetScan,
##'  miRDB and mirTarBase
##' @param DE.miRNA Character string or character vector of miRNA(s)
##' @param DE.target Character string or character vector of mRNA(s)
##' @param path Character string with the path to where the database files are stored
##' @return Returns a list with four elements: 1) list of  miRNAs each conataining 
##' a vector with all identified targets in DE.target; 2) dataframe with results from 
##' miRDB, 3) mirTarBase and4) TargetScan
##' @author Julia Perera Bel <jperera@imim.es>
##' @export
##' @import biomaRt
##' @import miRBaseConverter
##' @import gplots

miRNAGenes<-function (DE.miRNA, DE.target = NULL, path = "/bicoh/MARGenomics/annotationData/miRNAIntegration/"){
  library(org.Hs.eg.db)
  require(miRBaseConverter)
  require(gplots)
  
  # TargetScan: v8.0 Predicted.  Based on miRBase 22? n=228,052 interactions
  # miRDB: v6 (2019) Predicted. Based on miRBase 22. n=3,375,741 interactions
  # mirTarBase: v10 (2025) Validated. Based on miRBase 22. n=4,004,968 interactions
  
  #DE.miRNA: character string or character vector for the mature miRNA(s)
  #DE.target: character string or character vector for the target gene(s
  
  # load DBs
  message("Loading targetsacan...\n")
  targetscan_default=read.delim(file.path(path,"Predicted_Targets_Context_Scores.default_predictions_v8.0_HomoSapiens.txt"),header = T)
  message("Loading miRDB and mapping REFSEQ to SYMBOLS...\n")
  miRDB=read.delim(file.path(path,"miRDB_v6.0_prediction_result_HomoSapiens.txt"),header = F)
  miRDB$SYMBOL <- mapIds(org.Hs.eg.db, keys = miRDB$V2, column = "SYMBOL", keytype = "REFSEQ", multiVals = "first")
  message("Loading mirtarbase...\n")
  mirtarbase=read.delim(file.path(path,"miRTarBase_hsa_MTI_v10.csv"),header = T,quote = "",sep = ",")
  
  # Take as universe all target genes in the 3 databases
  if (is.null(DE.target)) {
    DE.target <- unique(c(mirtarbase$Target.Gene,targetscan_default$Gene.Symbol,miRDB$SYMBOL))
  }
  
  #######
  # miRDB (predicted)
  #######
  # Match with DB
  match1=miRDB[as.character(miRDB$V1) %in% DE.miRNA,]
  match1=match1[match1$SYMBOL %in% DE.target,]
  match1=match1[!(duplicated(match1[,c("V1","SYMBOL")])),]
  
  #######
  # mirtarbase (validated)
  #######
  # Match with DB
  match2=mirtarbase[mirtarbase$miRNA %in% DE.miRNA,]
  match2=match2[match2$Target.Gene %in% DE.target,]
  match2=match2[!(duplicated(match2[,c(2,4)])),]
  
  
  #######
  # targetscan
  #######
  # Match with DB "Default predictions (conserved sites of conserved miRNA families)"
  match3=targetscan_default[targetscan_default$miRNA %in% DE.miRNA,]
  match3=match3[match3$Gene.Symbol %in% DE.target,]

  v=venn(list("miRDB"=unique(match1$SYMBOL),
                  "miRtarBase"=unique(match2$Target.Gene),
                  "TargetScan"=unique(match3$Gene.Symbol)))
  plot(v)
  
  ##################
  #organize by miRNA
  gens.sel=list()
  for (i in DE.miRNA){
    gens.sel[[i]]<-unique(c(match1$SYMBOL[match1$V1==i],
                            as.character(match2$Target.Gene[match2$miRNA==i]),
                            as.character(match3$Gene.Symbol[match3$miRNA==i])))
    
  }
  return(list("miRNA-targets"=gens.sel,"miRDB"=match1,"miRTarBase"=match2,"TargetScan"=match3))
  
}
