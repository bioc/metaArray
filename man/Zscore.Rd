\name{Zscore}
\alias{Zscore}
\title{Meta-analysis of Microarray Data from Different Platforms}
\description{
  This function calculates Z-score for each matched gene across all datasets. 
  In each dataset, it performs local regression smoothing of mean vs variance.
  Z score is constructed by taking the ratio of weighted mean difference and 
  combined standard deviation according to Box and Tiao (1992).
}
\usage{
  Zscore(merged, pheno = NULL, permute = 0, verbose = TRUE)
}
\arguments{
  \item{merged}{\code{mergeExprSet} object that contains gene expression and class label with all
datasets. Class label should consist of two unique elements. If pheno is NULL, first columns
of phenoData of each \code{exprSet} is sought as class labels. If a vector of particular column number
in each data is specified, corresponding columns of phenoData will be considered for class labels. }
  \item{pheno}{A numeric vector specifying the location of class labels in phenoData slot of each 
\code{exprSet}, a unit of \code{mergeExprSet} representing one dataset. }
  \item{permute}{If permute is 0, weighted Z-score will be referenced to standard normal
distribution for two-sided p-value. Otherwise, columns of all datasets (each dataset separately) will be 
shuffled at random, from which a permutation distribution of Z-scores are formed and Z-scores are 
referenced to this distribution. }
  \item{verbose}{If verbose is TRUE, the progress of permutation will be reported. }
}
\value{
  A data.frame with matched genes, Z-scores and p-values will result.
}
\references{J.Wang et al, Bioinformatics 2004 Nov 22;20(17):3166-78  }
\author{Debashis Ghosh <ghoshd@umich.edu>, Hyungwon Choi <hwchoi@umich.edu> }
\examples{
 # Zscore(merged, pheno=NULL, permute=10000, verbose=FALSE)
}
\keyword{internal}
