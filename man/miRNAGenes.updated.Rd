\name{miRNAGenes.updated}
\alias{miRNAGenes.updated}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get microRNA-target Interactions.
}
\description{
  Retrieve predicted and validated miRNA-target interactions using TargetScan, miRDB and mirTarBase.
}
\usage{
miRNAGenes.updated(DE.miRNA, DE.target, path = "/bicoh/MARGenomics/annotationData/miRNAIntegration/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DE.miRNA}{
    Character string or character vector of miRNA(s)
}
  \item{DE.target}{
    Character string or character vector of mRNA(s)
}
  \item{path}{
    Character string with the path to where the database files are stored
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (DE.miRNA, DE.target, path = "/bicoh/MARGenomics/annotationData/miRNAIntegration/") 
{
    require(biomaRt)
    require(miRBaseConverter)
    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    DE.miRNA
    DE.target = as.character(DE.target)
    cat("Loading targetsacan...\n")
    targetscan_default = read.delim("/bicoh/MARGenomics/annotationData/miRNAIntegration/Predicted_Targets_Context_Scores.default_predictions_HomoSapiens.txt", 
        header = T)
    cat("Loading miRDB...\n")
    miRDB = read.delim("/bicoh/MARGenomics/annotationData/miRNAIntegration/miRDB_v6.0_prediction_result_HomoSapiens.txt", 
        header = F)
    cat("Loading mirtarbase...\n")
    mirtarbase = read.delim("/bicoh/MARGenomics/annotationData/miRNAIntegration/miRTarBase_MTI_HomoSapiens.txt", 
        header = T)
    v22 = miRNAVersionConvert(DE.miRNA, targetVersion = "v22", 
        exact = TRUE)
    refseq_mrna <- getBM(attributes = c("refseq_mrna", "hgnc_symbol"), 
        filters = "hgnc_symbol", values = as.character(DE.target), 
        mart = ensembl)
    match1 = miRDB[as.character(miRDB$V1) \%in\% v22$TargetName, 
        ]
    match1 = match1[match1$V2 \%in\% refseq_mrna$refseq_mrna, ]
    match1$V2 = refseq_mrna$hgnc_symbol[match(match1$V2, refseq_mrna$refseq_mrna)]
    match1 = match1[!(duplicated(match1[, 1:2])), ]
    match2 = mirtarbase[mirtarbase$miRNA \%in\% DE.miRNA, ]
    match2 = match2[match2$Target.Gene \%in\% DE.target, ]
    match2 = match2[!(duplicated(match2[, c(2, 4)])), ]
    match3 = targetscan_default[targetscan_default$miRNA \%in\% 
        DE.miRNA, ]
    match3 = match3[match3$Gene.Symbol \%in\% DE.target, ]
    venn(list(miRDB = match1$V2, miRtarBase = match2$Target.Gene, 
        TargetScan = match3$Gene.Symbol))
    gens.sel = unique(c(match1$V2, match2$Target.Gene, match3$Gene.ID))
    gens.sel = list()
    for (i in DE.miRNA) {
        gens.sel[[i]] <- unique(c(match1$V2[match1$V1 == i], 
            as.character(match2$Target.Gene[match2$miRNA == i]), 
            as.character(match3$Gene.Symbol[match3$miRNA == i])))
    }
    return(list(gens.sel, match1, match2, match3))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
