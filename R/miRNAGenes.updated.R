miRNAGenes.updated<-function(DE.miRNA,DE.target,path="/bicoh/MARGenomics/annotationData/miRNAIntegration/"){
  require(biomaRt)
  require(miRBaseConverter)
  require(gplots)
  ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
  # TargetScan: v7.2 2018 Predicted.  Based on miRBase 21
  # miRDB: v6 (2019) Predicted. Based on miRBase 22 !!! library(miRBaseConverter)
  # mirTarBase: v7 (2017) Validated. Based on miRBase 21
  
  #DE.miRNA: character string or character vector for the mature miRNA(s)
  #DE.tarfget: character string or character vector for the target gene(s
  DE.miRNA
  DE.target=as.character(DE.target)
  # load DBs
  cat("Loading targetsacan...\n")
  targetscan_default=read.delim(file.path(path,"Predicted_Targets_Context_Scores.default_predictions_HomoSapiens.txt"),header = T)
  cat("Loading miRDB...\n")
  miRDB=read.delim(file.path(path,"miRDB_v6.0_prediction_result_HomoSapiens.txt"),header = F)
  cat("Loading mirtarbase...\n")
  mirtarbase=read.delim(file.path(path,"miRTarBase_MTI_HomoSapiens.txt"),header = T)
  
  #######
  # miRDB
  #######
  
  # Convert IDs
  v22 = miRNAVersionConvert(DE.miRNA,targetVersion = "v22",exact = TRUE)
  refseq_mrna<-getBM(attributes=c("refseq_mrna","hgnc_symbol"),filters="hgnc_symbol",
                     values=as.character(DE.target),mart=ensembl) #from symbol to refseq (NM_)
  # Match with DB
  match1=miRDB[as.character(miRDB$V1) %in% v22$TargetName,]
  match1=match1[match1$V2 %in% refseq_mrna$refseq_mrna,]
  match1$V2=refseq_mrna$hgnc_symbol[match(match1$V2,refseq_mrna$refseq_mrna)]
  match1=match1[!(duplicated(match1[,1:2])),]
  
  #######
  # mirtarbase
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

  venn(list("miRDB"=match1$V2,"miRtarBase"=match2$Target.Gene,"TargetScan"=match3$Gene.Symbol))
  gens.sel=unique(c(match1$V2,match2$Target.Gene,match3$Gene.ID))
  
  
  ##################
  #organize by miRNA
  gens.sel=list()
  for (i in DE.miRNA){
    # venn(list("miRDB"=match1$V2[match1$V1==i],
    #           "miRtarBase"=as.character(match2$Target.Gene[match2$miRNA==i]),
    #           "TargetScan"=as.character(match3$Gene.Symbol[match3$miRNA==i])))
    # title(i)
    gens.sel[[i]]<-unique(c(match1$V2[match1$V1==i],
                            as.character(match2$Target.Gene[match2$miRNA==i]),
                            as.character(match3$Gene.Symbol[match3$miRNA==i])))
    
  }
  return(list(gens.sel,match1,match2,match3))
  
}