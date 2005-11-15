\name{intcor}
\alias{intcor}
\title{Integrative Correlation Analysis}
\description{
  This function calculates gene-specific reproducibility score
  based on Parmigiani et al. R implementation in MergeMaid (pairwise.cor) is a 
very efficient function in that it calculates the correlation matrix only 
once and collect appropriate elements for calculation of scores for each gene. 
However, in case there are more than thousands of common genes across datasets, 
the correlation matrix may overflow memory cells alloted to a session of R. 
Therefore, a replacement to the function that remedies the storage problem by brute 
force but fast computation in C is provided here. } 
\usage{
  intcor(merged)
}
\arguments{
  \item{merged}{\code{mergeExprSet} object that contains gene expression and class label with all
datasets.}
}
\value{
  \item{avg.cor}{A vector of gene-specific integrative correlation score}
  \item{pair.cor}{A matrix of correlations for each gene in every pair of two studies}
}
\references{ Clinical Cancer Research, Parmigiani et al. 10(9):2922-2927, 2004 }
\author{Debashis Ghosh <ghoshd@umich.edu>, Hyungwon Choi <hwchoi@umich.edu> }
\examples{
 #intcor(merged)
}
\keyword{internal}
