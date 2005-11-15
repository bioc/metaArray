\name{find.init}
\alias{find.init}
\title{Initialization of EM algorithm }
\description{
  This function is an automated initialization of 'z' in EM algorithm.
}
\usage{
  find.init(z, width = 1)
}
\arguments{
  \item{z}{ Unobserved probability of membership to uniform component of the mixture. }
  \item{width}{ Constant factor used when assigning 0/1 labels to samples. 
                Larger width will result in more samples initialized at 0. 
   } 
}
\value{
  \item{z}{A vector of 0/1, initial values of the EM algorithm }
}
\author{Debashis Ghosh <ghoshd@umich.edu>, Hyungwon Choi <hwchoi@umich.edu> }
\examples{
}
\keyword{internal}
