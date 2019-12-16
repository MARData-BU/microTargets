\name{miRNACorrel}
\alias{miRNACorrel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get microRNA-target Interactions and compute correlations.
}
\description{
    Main function to retrieve interactions between a set of miRNAs and a set of genes and computes correlations given expression data. The function calls miRNAGenes.updated function which uses TargetScan, miRDB and mirTarBase.
}
\usage{
miRNACorrel(DE.miRNA, DE.target, counts.miRNA, counts.mRNA, resultsDir = getwd(), name = "Summary")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DE.miRNA}{
    Character string or character vector of miRNA(s)
}
  \item{DE.target}{
    Character string or character vector of mRNA(s)
}
  \item{counts.miRNA}{
    Table of miRNA counts. First column contains miRNA ID's (miRBase v21)
}
  \item{counts.mRNA}{
    Table of mRNA counts with sample columns in teh same order as counts.miRNA. First column contains Symbol IDs.
}
  \item{cols.s}{
    Vector with colors corresponding to the conditions of each sample.
}
  \item{resultsDir}{
    Output directory.
}
  \item{name}{
    Name of the output files.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The function returns a dataframe with all miRNAs and the number of targets identified and the number of targets significantly negatively correlated for  each miRNA. It also generates:
  \item{Excel file }{One sheet for each miRNA, in which the correlations results for each target are shown (Rho and p.value of the spearman correlation)}
  \item{PDF }{One correlation plot for each miRNA-target pair}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Júlia Perera Bel <jperera@imim.es>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
