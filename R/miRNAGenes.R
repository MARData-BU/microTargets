miRNAGenes<-function(DE.miRNA,DE.target){
  require(openxlsx)
  require(multiMiR)
  require(miRNAtap)
  require(miRNAtap.db)
  require(biomaRt)
  require(SpidermiR)
  ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
  #DE.miRNA: character string or character vector for the mature miRNA(s)
  #DE.tarfget: character string or character vector for the target gene(s
  
  # load SpidermiR DBs to avoid doing it in every iteration
  validated_all<-SpidermiRdownload_miRNAvalidate(validated)
  circulating_all<-SpidermiRdownload_miRNAextra_cir(miRNAextra_cir)

  gens.sel=list()
  for (i in DE.miRNA){
    print(i)

    ############
    # Multimir #
    ############
    # Retrieve multimir interactions
    intersect <- get_multimir(org     = "hsa",
                              mirna   = i,
                              target  = DE.target,
                              table   = "all",
                              summary = TRUE,
                              predicted.cutoff.type = "p",
                              predicted.cutoff      = 10,
                              use.tibble = F)
    multimir=intersect@data
    #Minimum in 2 databases
    t=table(multimir$target_symbol,multimir$database)
    miRNA.genes.deg<-rownames(t)[rowSums(t!=0)>1]
    miRNA.genes.deg<-miRNA.genes.deg[miRNA.genes.deg!="NA"]
    if(length(miRNA.genes.deg)!=0){multimir=multimir[grep(paste(miRNA.genes.deg,collapse = "|"),multimir$target_symbol),]}
    if(length(miRNA.genes.deg)==0){multimir=multimir[0,]}

    ############
    # miRNAtap #
    ############

    mirnatap = getPredictedTargets(i, species ='hsa',method ='geom', min_src = 2) #min 2 sources
    miRNA.genes<-getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters="entrezgene_id",
                       values=rownames(mirnatap),mart=ensembl) #from enrez to symbols

    miRNA.genes.deg<-unique(intersect(miRNA.genes$hgnc_symbol,DE.target)) # targets also downreagulated
    entrez=miRNA.genes$entrezgene_id[miRNA.genes$hgnc_symbol %in% miRNA.genes.deg]
    mirnatap = mirnatap[as.character(entrez),]
    rownames(mirnatap)=miRNA.genes$hgnc_symbol[match(rownames(mirnatap),miRNA.genes$entrezgene_id)]
    colnames(mirnatap)[1:5]=c('pictar','diana','targetscan','miranda','mirdb')

    #############
    # SpidermiR #
    #############

    pred=SpidermiRdownload_miRNAprediction(mirna_list=i)
    pred=pred[grep(paste(miRNA.genes.deg,collapse = "|"),pred$V2),]
    pred$type=rep("predicted",nrow(pred))

    val=validated_all[validated_all$V1==i,]
    val=val[grep(paste(miRNA.genes.deg,collapse = "|"),val$V2),]
    val$type=rep("validated",nrow(val))

    #We are not using this information....
    circulating_all=circulating_all[circulating_all$organism=="Homo sapiens",]
    circ=circ[circ$miRBase_Last_Version==i,]

    spidermir=rbind(val,pred)
    colnames(spidermir)=c("miRNA","target","type")

    ##################
    #select those at least retrieved by 2 packages
    tots.g<-unique(c(multimir$target_symbol,rownames(mirnatap),as.character(spidermir$target)))
    lng=length(tots.g)
    gens.sel.symbol<- ""
    if (lng>0) {
      tots.g.n<-array(0,dim=lng)
      for (j in 1:lng){
        if (tots.g[j] %in% multimir$target_symbol) tots.g.n[j]=tots.g.n[j]+1
        if (tots.g[j] %in% rownames(mirnatap)) tots.g.n[j]=tots.g.n[j]+1
        if (tots.g[j] %in% as.character(spidermir$target)) tots.g.n[j]=tots.g.n[j]+1
      }
      #select those at least retrieved by 2 packages
      gens.sel[[i]]<-tots.g[which(tots.g.n>2)]
    }
  }

  return(gens.sel)
}
