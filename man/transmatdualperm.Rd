\name{transmatdualperm}
\alias{transmatdualperm}
%\alias{transmatdualperm2}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{transmatdualperm}
\description{
Works out the transform matrix for a particular prediction scheme and neighbourhood structure.
}
\usage{
transmatdualperm(x, f, Pred = AdaptNeigh, neigh = 1, int = TRUE, clo = 
TRUE, keep = 2,perm = 
sample(1:length(x),(length(x)-keep),FALSE),varonly=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of grid values. Can be of any length, not necessarily equally spaced.}
  \item{f}{A vector of function values corresponding to \env{x}. Must be of the same length as \env{x}.}
  \item{Pred}{The type of regression to be performed. Possible options are \command{\link{LinearPred}}, \command{\link{QuadPred}}, \command{\link{CubicPred}}, \command{\link{AdaptPred}} and \command{\link{AdaptNeigh}}.}
  \item{neigh}{The number of neighbours over which the regression is performed at each step. If \env{clo} is false, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{int}{Indicates whether or not the regression curve includes an intercept.}
  \item{clo}{Refers to the configuration of the chosen neighbours. If \env{clo} is false, the neighbours will be chosen symmetrically around the removed point. Otherwise, the closest neighbours will be chosen.}
  \item{keep}{The number of scaling coefficients to be kept in the final representation of the initial signal. This must be at least two.}
  \item{perm}{Vector of length (length(\code{x})-\code{keep}). It gives the trajectory for the modified lifting algorithm to follow, i.e. it gives 
  the order of point removal.}
  \item{varonly}{A boolean variable indicating whether only the 
  coefficient variances should be returned, i.e. just the diagonal of 
  the transform matrix \code{Wnew}.}
}
\details{
 The function uses \command{Amatdual} to form the refinement matrices \eqn{A_j}, from which the augmented matrices \eqn{T_j} are constructed. This process is iterated, to find the transform matrix (the top level augmented matrix). The rows and columns of this matrix are then reordered to be in the order of \code{out$coeff}, i.e. so that the columns correspond to \eqn{f_1 \dots f_n}.
}
\value{
    \item{out}{the output from the forward transform.}
  \item{Wnew}{the matrix associated to the wavelet transform.}
  \item{x}{the original gridpoint vector. }
}

\author{Matt Nunes (\email{m.nunes@lancs.ac.uk}), Marina Knight }

\note{ This function has been left in the package for completeness.  However, the transform matrix is (optionally) computed within the forward lifting
transform function \code{\link{fwtnp}}.
}

\seealso{
\code{\link{fwtnp}}
}

\examples{
x1<-runif(10)
y1<-make.signal2("doppler",x=x1)
#
vec<-sample(1:10,8,FALSE)

a<-transmatdualperm(x1,y1,AdaptNeigh,2,TRUE,TRUE,2,perm=vec)
#
a$Wnew
#
#the transform matrix for this adaptive lifting scheme 

}
\keyword{array}% at least one, from doc/KEYWORDS


