\name{poe.em}
\alias{poe.em}
\title{Probability of Expression from mixture distribution for multiple genes. }
\description{
  This function applies fit.em function to all rows of a gene expression data matrix. 
}
\usage{
  poe.em(mat, cl = NULL, threshold=1e-05, every = 100)
}
\arguments{
  \item{mat}{Gene expression data matrix }
  \item{cl}{A vector of 0s and 1s. Use 1 for normal phenotype and 0 for
    non-normal phenotype. If left blank, all samples will be labeled as
    normal phenotype. Normal component of mixture is estimated using
    samples with normal phenotype only. POE for samples with non-normal
    phenotype will be calculated after EM algorithm finishes ML estimation.}
  \item{threshold}{Criterion of convergence in likelihood }
  \item{every}{Progress of estimation is reported at every integer mode of the value 'every' }
}
\value{
  A data matrix of transformed expression will result.
}
\author{Debashis Ghosh <ghoshd@umich.edu>, Hyungwon Choi <hwchoi@umich.edu> }
\examples{}
\keyword{internal}
