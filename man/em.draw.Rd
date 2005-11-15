\name{em.draw}
\alias{em.draw}
\title{Plot of transformed expression produced by EM algorithm }
\description{
  Given a numeric vector, a plot of four panels is drawn:
    1) fitted mixture distribution
    2) transformed expression against original expression
    3) histogram of original expression 
    4) progression of log-likelihood during the fit
}
\usage{
  em.draw(vec, cl=NULL, threshold=0.0001)
  #em.draw(vec, cl=1-metastasis, threshold=0.0001)
}
\arguments{
  \item{vec}{A numeric vector, especially a particular row of expression matrix }
  \item{cl}{A vector of 0s and 1s. Use 1 for normal phenotype and 0 for
    non-normal phenotype. If left blank, all samples will be labeled as
    normal phenotype. Normal component of mixture is estimated using
    samples with normal phenotype only. POE for samples with non-normal
    phenotype will be calculated after EM algorithm finishes ML estimation. }
  \item{threshold}{Crieterion for convergence in likelihood. }
}
\value{
  A plot of four panels will appear upon the call.
}
\author{Debashis Ghosh <ghoshd@umich.edu> , Hyungwon Choi <hwchoi@umich.edu> }
\examples{

}
\keyword{internal}
