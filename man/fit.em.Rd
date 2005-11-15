\name{fit.em} 
\alias{fit.em} 
\title{Probability of expresison from mixture distribution for a single gene. } 
\description{
  This core function fits two-component normal-uniform mixture
  distribution, and extracts probability of over/under expression for all samples in all genes.  
}
\usage{
  fit.em(x, cl, threshold=1e-05)
}
\arguments{
  \item{x}{A numeric vector, especially expression values for a particular gene. }
  \item{cl}{A vector of 0s and 1s. Use 1 for normal phenotype and 0 for
    non-normal phenotype. If left blank, all samples will be labeled as
    normal phenotype. Normal component of mixture is estimated using
    samples with normal phenotype only. POE for samples with non-normal
    phenotype will be calculated after EM algorithm finishes ML estimation.}
  \item{threshold}{Criterion for convergence in likelihood }
}
\value{
  \item{expr }{Estimated POE }
  \item{a }{Minimum (adjusted) of Raw Expression }
  \item{b }{Maximum (adjusted) of Raw Expression }
  \item{sigmasq }{Estimated variance of normal component }
  \item{mu }{Estimated mean of normal component }
  \item{pi }{Probability that the gene is over/under expressed on
    average across the samples}
  \item{conv }{Convergence status of EM }
  \item{dec }{Detection of decreasing likelihood }
  \item{loc }{Signs of local maxima }
  \item{lik.rec }{Trajectory of likelihood during EM }
  \item{abnormal }{Overall sign of failure in EM by any of the three criteria above }
}
\author{Debashis Ghosh <ghoshd@umich.edu>, Hyungwon Choi <hwchoi@umich.edu> }
\examples{

}
\keyword{internal}
