
\name{poe.mcmc}
\alias{poe.mcmc}
\title{Probability of Expression (POE)}
\description{
  Differential expression using latent categories of up/down regulation. 
  Three component normal-uniform mixture model under Bayesian hierarchical analysis.
  This is a C implementation of poe.fit function of POE package (MCMC). R portion of source code
  was directly adapted from POE; poe.one.iteration function was re-written in C.
  Some of the optional arguments available in poe.fit are suppressed here, and you cannot save the chain of
  samples drawn for numerical integration.
}
\usage{
  poe.mcmc(AA, NN = NULL, id = NULL, M = 2000, kap.min=3.0,
      logdata=FALSE, stepsize=0.5, 
      centersample = FALSE, centergene = FALSE, generatestarts = TRUE, start.method = 1,
      startobject = R0, collapse.to.two = FALSE, burnin=200,
      collapse.window=50, converge.threshold=0.01,
      PR = list(alpha.mm = 0, alpha.sd = 100, mu.mm = 0, mu.sd = 100,
          pipos.mm = 0, pipos.sd = 100, pineg.mm = 0, pineg.sd = 100,
          kap.pri.rate = 1, tausqinv.aa = 1, tausqinv.bb = 0.1))
  # poe.mcmc(AA = chen, NN = 1 - chen.spl$metastasis, M=2000)
}
\arguments{
  \item{AA}{Matrix or \code{exprs} from an \code{ExpressionSet} object.}
  \item{NN}{Phenotypic label of arrays. Samples from the reference group (normal phenotype) 
            must be labeled as 1. For instance, if you wish to identify differentially expressed
            genes between cancer and normal samples, mark normal as 1 and cancer as 0. Then the 
            high probability hits indicate over or under-expressed genes in specific samples 
            with respect to the gene expression distribution of the normal samples.
            Alternatively, you may leave this blank if you wish to get the probability of high and low
            expression from the mixture distribution using all samples not relying on the known sample information.}
  \item{M}{Number of MCMC iterations after burn-in.}
}
\value{
  \item{poe }{Probability of over/under expression. Transformed gene expression on a fixed 
scale of [-1,1].}
  \item{Other values }{Posterior median estimates of parameters. Please refer to POE package for details.}
}
\references{G. Parmigiani et al, JRSS, 64:717-736, 2002 or URL: http://astor.som.jhmi.edu/poe/}
\author{Debashis Ghosh <ghoshd@psu.edu>, Hyungwon Choi <hyung_won_choi@nuhs.edu.sg>}
\examples{
  # poe.mat <- poe.mcmc(AA=exprmat, NN=clvec, M=10000)
  # One can also provide different hyperparameter values. 
}
\keyword{internal}
