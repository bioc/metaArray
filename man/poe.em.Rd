\name{poe.em}
\alias{poe.em}
\title{Probability of Expression from mixture distribution for multiple genes. }
\description{
  This function applies fit.em function to all rows of a gene expression data matrix. 
}
\usage{
  poe.em(mat, cl, threshold=1e-05, every = 100)
}
\arguments{
  \item{mat}{Gene expression data matrix }
  \item{cl}{A vector of 0s and 1s. Use 1 for normal phenotype and 0 for
    non-normal phenotype. If all samples are of unknown phenotype,
    give vectors of 0. This is identical to the poe.mcmc function. }
  \item{threshold}{Criterion of convergence in likelihood }
  \item{every}{Progress of estimation is reported at every integer mode of the value 'every' }
}
\value{
  A data matrix of transformed expression will result.
}
\author{Debashis Ghosh <ghoshd@psu.edu>, Hyungwon Choi <hyung_won_choi@nuhs.edu.sg> }
\examples{}
\keyword{internal}
