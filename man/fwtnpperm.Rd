\name{fwtnpperm}
\alias{fwtnpperm}
\alias{fwtnpperm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fwtnpperm}
\description{
Performs the lifting transform on a signal with grid \env{input} and corresponding function values \env{f}. There is a unique correspondence between the grid values and the function values. Can also cope with length vector input instead of gridpoint vector input.
}
\usage{
fwtnpperm (input, f, LocalPred = LinearPred, neighbours = 1,
intercept = TRUE,closest = FALSE, nkeep = 2, initboundhandl = "reflect", mod =
    sample(1:length(input), (length(input) - nkeep), FALSE),
    do.W = FALSE, varonly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{A vector of grid values. Can be of any length, not necessarily equally spaced.}
  \item{f}{A vector of function values corresponding to \env{input}. Must be of the same length as \env{input}.}
  \item{LocalPred}{The type of regression to be performed. Possible options are \code{LinearPred}, \code{QuadPred}, \code{CubicPred}, \code{AdaptPred} and \code{AdaptNeigh}.}
  \item{neighbours}{The number of neighbours over which the regression is performed at each step. If closest is false, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{intercept}{Indicates whether or not the regression curve includes an intercept.}
  \item{closest}{Refers to the configuration of the chosen neighbours. If \env{closest} is false, the neighbours will be chosen symmetrically around the removed point. Otherwise, the closest neighbours will be chosen.}
  \item{nkeep}{The number of scaling coefficients to be kept in the final representation of the initial signal. This must be at least two. }
  \item{initboundhandl}{variable specifying how to handle the boundary at the start of the transform.  Possible values are \code{"reflect"} - the intervals corresponding to the first and last datapoints are taken to have the respective grid values as midpoints; and \code{"stop"} - the first and last intervals have the first and last grid values (respectively) as outer endpoints.}
  \item{mod}{Vector of length (length(\code{x})-\code{keep}). It gives 
  the trajectory for the modified lifting algorithm to follow, i.e. it gives 
  the order of point removal.}
  \item{do.W}{A boolean indicating whether the transform matrix should 
  be computed and returned.}
  \item{varonly}{A boolean indicating whether only the coefficient 
  variances should be returned (if \code{do.W=TRUE}).}
}
\details{
Given \eqn{n} points on a line, \env{input}, each with a corresponding env{f} value this algorithm computes a lifting transform of the (\env{input},\env{f}) data. If lengths are inputted (\env{inputtype="lengths"}), then the gridpoints are taken to be the left endpoints of the intervals defined by the lengths inputted. 
Step One. Order the grid values so that corresponding intervals can be constructed. 

Step Two. Compute "integrals" for each point. For each point its integral is the length of the interval associated to the gridpoint. 

Step Three. Identify the point to remove as that with the smallest integral. Generally, we remove points in order of smallest to largest integral. The integrals of neighbours of removed points change at each step. 

Step Four(a). The neighbours of the removed point are identified using the specified neighbour configuration. The value of \env{f} at the removed point is predicted using the specified regression curve over the neighbours, unless an adaptive procedure is chosen. In this case, the algorithm adjusts itself. The difference between the removed point's \env{f} value and the prediction is computed: this is the wavelet coefficient for the removed point. The difference replaces the function value in the vector coeff at the removed point's location. In this way wavelet coefficients gradually overwrite (scaling) function values in \env{coeff}. 

Step Four(b). The integrals and the scaling function values (other \env{coeff} values) of neighbours of the removed point are updated. The values of the rest of the scaling coefficients are unaffected. 

Step Five. Return to step 3 but in the identification of a point to remove the updated integrals are used. 

The algorithm continues until as many points as desired are removed.  If 
\code{do.W=TRUE}, the predict and update lifting steps are used to 
propogate coefficient contributions to the transform matrix \code{W}.  
If \code{varonly=TRUE}, only the (detail and scaling) coefficient 
variances are returned.  After each lifting step, the 
coefficient variance is computed and the transform matrix 
row corresponding to the lifted coefficient is deleted for the 
next stage (minimal storage efficiency).  The transform matrix is not 
returned (i.e. \code{W=NULL}). 
}
\value{
\item{X}{data vector of the grid used in the transform.}
\item{coeff}{vector of detail and scaling coefficients in the wavelet decomposition of the signal.}
\item{origlengths}{vector of initial interval lengths corresponding to the gridpoints.}
\item{lengths}{vector of (updated) interval lengths at the end of the transform. This is of length \env{nkeep}.}
\item{lengthsremove}{vector of interval lengths corresponding to the points removed during the transform (in \env{removelist}).}
\item{pointsin}{indices into \env{X} of the scaling coefficients in the wavelet decomposition. These are the indices of the \env{X} values which remain after all points in \env{removelist} have been predicted and removed. This has length \env{nkeep}.}
\item{removelist}{a vector of indices into \env{X} of the lifted coefficients during the transform (in the order of removal).}
\item{neighbrs}{a list of indices into \env{X}. Each list entry gives the indices of the neighbours of the removed point used at that particular step of the transform.}
\item{neighbours}{the user-specified number of neighbours used in the prediction step of the transform.}
\item{gamlist}{a list of all the prediction weights used at each step of the transform.}
\item{alphalist}{a list of the update coefficients used in the update step of the decomposition.}
\item{schemehist}{a vector of character strings indicating the type of regression used at each step of the transform.}
\item{interhist}{a boolean vector indicating whether or not an intercept was used in the regression curve at each step.}
\item{clolist}{a boolean vector showing whether or not the neighbours were symmetrical (FALSE) about the removed point during the transform. This is NULL except when \code{LocalPred=AdaptNeigh}.}
}

\author{Matt Nunes (\email{m.nunes@lancs.ac.uk}), Marina.Knight}

\seealso{\code{\link{AdaptNeigh}}, \code{\link{AdaptPred}}, \code{\link{CubicPred}}, \code{\link{fwtnpmp}}, \code{\link{invtnp}}, \code{\link{LinearPred}}, \code{\link{QuadPred}}     }
\examples{
#
# Generate some one-dimensional data: 100 observations.
#
input <- runif(100)
f <- input^2 - 3*input
#
# Compute fwtnp function on this data
#
vec<-sample(1:100,98,FALSE)

out <- fwtnpperm(input,f,LocalPred=AdaptPred,neighbours=2,closest=TRUE,mod=vec)
#
# That's it.
#
}
\keyword{methods}% at least one, from doc/KEYWORDS

