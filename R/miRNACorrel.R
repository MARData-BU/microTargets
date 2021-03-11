##' Get microRNA-target Interactions and compute correlations
##'
##' Main function to retrieve interactions between a set of miRNAs and a 
##' set of genes and computes correlations given expression data. The function 
##' calls miRNAGenes.updated function which uses TargetScan, miRDB and mirTarBase.
##' @param DE.miRNA Character string or character vector of miRNA(s)
##' @param DE.target Character string or character vector of mRNA(s)
##' @param counts.miRNA Table of miRNA counts. First column contains miRNA ID's (miRBase v21)
##' @param counts.mRNA Table of mRNA counts with sample columns in teh same 
##' order as counts.miRNA. First column contains Symbol IDs
##' @param cols.s Vector with colors corresponding to the conditions of each sample
##' @param resultsDir Output directory (default working dir)
##' @param name Name of the output files (default "Summary")
##' @param path Character string with the path to where the database files are stored
##' @return The function returns a dataframe with all miRNAs and the number of 
##' targets identified and the number of targets significantly negatively 
##' correlated for  each miRNA. It also generates: a) excel file with one sheet 
##' per miRNA, in which the correlations results for each target are shown 
##' (Rho and p.value of the spearman correlation); b) pdf with correlation plots 
##' for each miRNA-target pair
##' @author Julia Perera Bel <jperera@imim.es>
##' @export
##' @import openxlsx
miRNACorrel<-function(DE.miRNA,DE.target, counts.miRNA, counts.mRNA,cols.s=NULL,resultsDir=getwd(),name="Summary"){
  require(openxlsx)

  #DE.miRNA: character string or character vector for the mature miRNA(s)
  #DE.tarfget: character string or character vector for the target gene(s
  #counts.miRNA: table of counts with Symbol in first column and columns of samples of interest in correct order
  #counts.mRNA: table of counts with miRNAs in rownames and only columns of samples of interest in correct order

  # For plotting if cols.s is not provided
  if (length(cols.s)==0){
    cols.s<-rep("green",ncol(counts.miRNA)-1)
  }

  #Get targets of miRNAs using 3 databases. Intersect reults with DE genes
  #gens.sel=miRNAGenes(DE.miRNA,DE.target)
  gens.sel=miRNAGenes.updated(DE.miRNA,DE.target)[[1]]

  miRNAs<-unique(names(gens.sel))
  wb <- createWorkbook()
  #wb <- loadWorkbook()
  cor.summary=data.frame(matrix(NA,nrow = length(miRNAs)))
  rownames(cor.summary)=miRNAs
  for (i in miRNAs){
    print(i)
    #miRNA.genes.deg<- unique(multimir[multimir$mature_mirna_id==i,"target_symbol"])# targets also downreagulated
    y=as.vector(counts.miRNA[i,]) # expr. data miRNA
    y=as.numeric(y)
    #now correlations
    miRNA.genes.deg=gens.sel[[i]]
    lng<-length(miRNA.genes.deg)
    cat(lng,"target genes")
    cor.rho<-array(NA,lng)
    names(cor.rho)=miRNA.genes.deg
    cor.pval<-array(NA,lng)
    names(cor.pval)=miRNA.genes.deg
    pdf(file.path(paste(resultsDir,paste(i,".corr.mRNA.miRNA.pdf"),sep="/")))
    for (j in miRNA.genes.deg){
      #print(j)
      mRNA<-j
      x=counts.mRNA[counts.mRNA$Symbol==mRNA,-1] #expr. data mRNA
      if(nrow(x)>1){ # In case two probes match to the same gene, we take gene with higher CV
        coeff.var=function(x){sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)}
        cv=apply(x,1,coeff.var)
        x=x[names(which(cv==max(cv))),]
      }
      x=as.numeric(x)
      cor<-cor.test(x,y, method = "spearman",exact=FALSE,alternative = "less") #busquem nomes cor. negatives
      cor.pval[j]<-cor$p.value
      cor.rho[j]<-cor$estimate
      plot(x,y,main=mRNA,xlab="log2RMA expression",ylab="log2miRMA expression",type="p",cex=0.8,col=cols.s,pch=19)
      legend("topright",legend = c("Control","Treatment"),fill = unique(cols),
             title = paste0("Rho: ",round(cor$estimate,2)," (p=",round(cor$p.value,2),")"),bty = "n")
      fit <- lm(y ~ x)
      abline(fit, col="chartreuse3")
    }
    dev.off()  #close pdf file
    cor.table<-data.frame("miRNA ID"=rep(i,lng),miRNA.genes.deg,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
    cor.summary[i,"num.targets"]=lng
    cor.summary[i,"rho.less.-0.6 & p<0.05"]=nrow(cor.table[cor.table$Rho<(-0.6) & cor.table$pval<0.05,])
    #cor.summary[i,"%.pval.less.0.05"]=round(100*sum(cor.table$pval<0.05)/lng,2)
    addWorksheet(wb, sheetName = i )
    writeData(wb,sheet = i,cor.table,startRow = 1, startCol = 1)
  }
  cor.summary=cor.summary[,-1]
  saveWorkbook(wb, paste(resultsDir,"/",name,"_corr.mRNA.miRNA.xlsx",sep=""), overwrite = T)
  return(cor.summary)

}


miRNACorrel.multimir<-function(DE.miRNA,DE.target, counts.miRNA, counts.mRNA,resultsDir=getwd(),name="Summary"){
  require(openxlsx)
  require(multiMiR)
  #DE.miRNA: character string or character vector for the mature miRNA(s)
  #DE.tarfget: character string or character vector for the target gene(s
  #counts.miRNA: table of counts with Symbol in first column and columns of samples of interest in correct order
  #counts.mRNA: table of counts with miRNAs in rownames and only columns of samples of interest in correct order

  # Retrieve multimir interactions
  intersect <- get_multimir(org     = "hsa",
                            mirna   = DE.miRNA,
                            target  = DE.target,
                            table   = "all",
                            summary = TRUE,
                            predicted.cutoff.type = "p",
                            predicted.cutoff      = 10,
                            use.tibble = F)

  #multimir=intersect@data[grep("mirecords|mirtarbase|tarbase|targetscan|mirdb",intersect@data$database),]
  multimir=intersect@data

  miRNAs<-unique(multimir$mature_mirna_id)
  wb <- createWorkbook()
  #wb <- loadWorkbook()
  cor.summary=data.frame(matrix(NA,nrow = length(miRNAs)))
  rownames(cor.summary)=miRNAs
  for (i in miRNAs){
    print(i)
    multimir.i=multimir[multimir$mature_mirna_id==i,]
    #Minimum in 3 databases
    t=table(multimir.i$target_symbol,multimir.i$database)
    miRNA.genes.deg<-rownames(t)[rowSums(t!=0)>1]
    miRNA.genes.deg<-miRNA.genes.deg[miRNA.genes.deg!="NA"]
    #miRNA.genes.deg<- unique(multimir[multimir$mature_mirna_id==i,"target_symbol"])# targets also downreagulated
    y=as.vector(counts.miRNA[i,]) # expr. data miRNA
    y=as.numeric(y)
    #now correlations
    lng<-length(miRNA.genes.deg)
    cat(lng,"target genes")
    cor.rho<-array(NA,lng)
    names(cor.rho)=miRNA.genes.deg
    cor.pval<-array(NA,lng)
    names(cor.pval)=miRNA.genes.deg
    pdf(file.path(paste(resultsDir,paste(i,".corr.mRNA.miRNA.pdf"),sep="/")))
    for (j in miRNA.genes.deg){
      #print(j)
      mRNA<-j
      x=counts.mRNA[counts.mRNA$Symbol==mRNA,-1] #expr. data mRNA
      if(nrow(x)>1){ # In case two probes match to the same gene, we take gene with higher CV
        coeff.var=function(x){sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)}
        cv=apply(x,1,coeff.var)
        x=x[names(which(cv==max(cv))),]
      }
      x=as.numeric(x)
      cor<-cor.test(x,y, method = "spearman",exact=FALSE,alternative = "less") #busquem nomes cor. negatives
      cor.pval[j]<-cor$p.value
      cor.rho[j]<-cor$estimate
      plot(x,y,main=mRNA,xlab="log2RMA expression",ylab="log2miRMA expression",type="p",cex=0.8,col=cols)
      legend("topright",legend = c("Ctrl","Tto"),fill = unique(cols))
      fit <- lm(y ~ x)
      abline(fit, col="chartreuse3")
    }
    dev.off()  #close pdf file
    cor.table<-data.frame("miRNA ID"=rep(i,lng),miRNA.genes.deg,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
    cor.summary[i,"num.targets"]=lng
    cor.summary[i,"%.rho.less.-0.6"]=round(100*sum(cor.table$Rho<(-0.6))/lng,2)
    cor.summary[i,"%.pval.less.0.05"]=round(100*sum(cor.table$pval<0.05)/lng,2)
    addWorksheet(wb, sheetName = i )
    writeData(wb,sheet = i,cor.table,startRow = 1, startCol = 1)
  }
  cor.summary=cor.summary[,-1]
  saveWorkbook(wb, paste(resultsDir,"/",name,"_corr.mRNA.miRNA.xlsx",sep=""), overwrite = T)
  return(list(multimir,cor.summary))

}


miRNACorrel.mirnatap<-function(DE.miRNA,DE.target, counts.miRNA, counts.mRNA,resultsDir=getwd(),name="Summary"){
  require(openxlsx)
  require(miRNAtap)
  require(miRNAtap.db)
  require(biomaRt)
  ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")

  #DE.miRNA: character string or character vector for the mature miRNA(s)
  #DE.tarfget: character string or character vector for the target gene(s
  #counts.miRNA: table of counts with Symbol in first column and columns of samples of interest in correct order
  #counts.mRNA: table of counts with miRNAs in rownames and only columns of samples of interest in correct order

  miRNAs<-unique(DE.miRNA)

  wb <- createWorkbook()
  #wb <- loadWorkbook()
  cor.summary=data.frame(matrix(NA,nrow = length(miRNAs)))
  rownames(cor.summary)=miRNAs
  for (i in miRNAs){
    print(i)

    #Get target genes using miRNAtap
    predictions = getPredictedTargets(i, species ='hsa',method ='geom', min_src = 2) #min 3 sources
    miRNA.genes<-getBM(attributes="hgnc_symbol",filters="entrezgene_id",
                       values=rownames(predictions),mart=ensembl)$hgnc_symbol #from enrez to symbols
    miRNA.genes.deg<-unique(intersect(miRNA.genes,DE.target)) # targets also downreagulated

    y=as.vector(counts.miRNA[i,]) # expr. data miRNA
    y=as.numeric(y)
    #now correlations
    lng<-length(miRNA.genes.deg)
    cat(lng,"target genes")
    cor.rho<-array(NA,lng)
    names(cor.rho)=miRNA.genes.deg
    cor.pval<-array(NA,lng)
    names(cor.pval)=miRNA.genes.deg
    pdf(file.path(paste(resultsDir,paste(i,".corr.mRNA.miRNA.pdf"),sep="/")))
    for (j in miRNA.genes.deg){
      #print(j)
      mRNA<-j
      x=counts.mRNA[counts.mRNA$Symbol==mRNA,-1] #expr. data mRNA
      if(nrow(x)>1){ # In case two probes match to the same gene, we take gene with higher CV
        coeff.var=function(x){sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)}
        cv=apply(x,1,coeff.var)
        x=x[names(which(cv==max(cv))),]
      }
      x=as.numeric(x)
      cor<-cor.test(x,y, method = "spearman",exact=FALSE,alternative = "less") #busquem nomes cor. negatives
      cor.pval[j]<-cor$p.value
      cor.rho[j]<-cor$estimate
      plot(x,y,main=mRNA,xlab="log2RMA expression",ylab="log2miRMA expression",type="p",cex=0.8,col=cols)
      legend("topright",legend = c("Ctrl","Tto"),fill = unique(cols))
      fit <- lm(y ~ x)
      abline(fit, col="chartreuse3")
    }
    dev.off()  #close pdf file
    cor.table<-data.frame("miRNA ID"=rep(i,lng),miRNA.genes.deg,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
    cor.summary[i,"num.targets"]=lng
    cor.summary[i,"%.rho.less.-0.6"]=round(100*sum(cor.table$Rho<(-0.6))/lng,2)
    cor.summary[i,"%.pval.less.0.05"]=round(100*sum(cor.table$pval<0.05)/lng,2)
    addWorksheet(wb, sheetName = i )
    writeData(wb,sheet = i,cor.table,startRow = 1, startCol = 1)
  }
  cor.summary=cor.summary[,-1]
  saveWorkbook(wb, paste(resultsDir,"/",name,"_corr.mRNA.miRNA.xlsx",sep=""), overwrite = T)
  return(cor.summary)

}
