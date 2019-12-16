\name{miRNAGenes}
\alias{miRNAGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get microRNA-target Interactions.
}
\description{
  Retrieve predicted and validated miRNA-target interactions using Multimir, miRNAtap and SpideRmiR R packages.
}
\usage{
  miRNAGenes(DE.miRNA, DE.target)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DE.miRNA}{
    Character string or character vector of miRNA(s)
  }
  \item{DE.target}{
    Character string or character vector of mRNA(s)
  }
}
\details{
  Requires several packages: multiMiR, miRNAtap, miRNAtap.db and SpidermiR. Not recomended to use this function as the database versions used in these packages are not up to date.
}
\value{
  Returns a list of miRNAs, each conataining a vector with all identified targets.

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
