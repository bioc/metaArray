\name{fit.em} 
\alias{fit.em} 
\title{Probability of expresison from mixture distribution for a single gene. } 
\description{
  This core function fits two-component normal-uniform mixture
  distribution, and extracts probability of over/under expression for all samples in all genes.  
}
\usage{
  fit.em(x, cl, threshold=1e-06)
}
\arguments{
  \item{x}{A numeric vector, especially expression values for a particular gene. }
  \item{cl}{A vector of 0s and 1s. Use 1 for normal phenotype and 0 for
    non-normal phenotype. Note that this is the opposite of POE MCMC.
    If all samples are of unknown phenotype or of the same one, give
    vector of zeros. When class information is provided, conditional estimation of the mixture 
    is applied. }
  \item{threshold}{Criterion for convergence in likelihood }
}
\value{
  \item{expr }{Estimated POE }
  \item{a }{Minimum (adjusted) of Raw Expression }
  \item{b }{Maximum (adjusted) of Raw Expression }
  \item{sigmasq }{Estimated variance of normal component }
  \item{mu }{Estimated mean of normal component }
  \item{Pi }{Probability that the gene is over/under expressed on
    average across the samples}
  \item{lik.rec }{Trajectory of likelihood during EM }
}
\author{Debashis Ghosh <ghoshd@umich.edu>, Hyungwon Choi <hwchoi@umich.edu> }
\examples{
  
}
\keyword{internal}
