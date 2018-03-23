\name{nlt}
\alias{nlt}
\alias{nlt2}
\title{Denoise a signal using a nondecimated lifting transform}
\description{
  Starting with a noise-contaminated signal, we decompose it using a 'nondecimated' lifting algorithm (i.e. by applying the modified lifting transform following several random paths), shrink all the obtained detail coefficients and invert each transform to give an estimated signal. The average of all these estimates is the final proposal for estimating the true (unknown) signal.
}
\usage{
nlt(x, f, J, Pred = AdaptPred, neighbours = 1, closest = FALSE, intercept = TRUE, nkeep = 2, trule = "median",verbose = TRUE,do.orig = FALSE, returnall = FALSE)
}
\arguments{
  \item{x}{Vector of any length (possibly irregularly spaced) that gives the grid locations where the signal is observed.}
  \item{f}{Vector of the same length as \code{x} that gives the signal values corresponding to the \code{x}-locations.}
  \item{J}{Number of trajectories to be used by the nondecimated lifting algorithm.}
  \item{Pred}{The type of regression to be used in the prediction step of the modified lifting algorithm. Choices are linear, quadratic or cubic (respectively, \code{LinearPred}, \code{QuadPred} or \code{CubicPred}), or two adaptive procedures which automatically choose the degree used in regression, (\code{AdaptPred} or \code{AdaptNeigh}).}
  \item{neighbours}{Number of neighbours to be used for defining the neighbourhood of each point that has to be removed. If (\code{closest}=FALSE), then this gives the number of neighbours to be used on each side of the removed point.}
  \item{closest}{If (\code{closest}=TRUE) or (\code{closest}=FALSE), then at each step the neighbours are in closest, respectively symmetrical configuration.}
  \item{intercept}{Specifies whether (\code{intercept}=TRUE) or not (\code{intercept}=FALSE) an intercept is to be used in the regression curve. For \code{Pred}=\code{AdaptPred} or \code{AdaptNeigh}, the algorithm automatically makes this choice.}
  \item{nkeep}{Number of scaling points we want at the end of the application of the transform. The usual choice is \code{nkeep}=2.}
  \item{trule}{The type of Bayesian shrinkage technique, with possible choices posterior median (\code{"median"}) or posterior mean (\code{"mean"}).}
  \item{verbose}{A boolean indicating whether extra information should be printed.}
  \item{do.orig}{A boolean indicating whether the original \code{adlift} algorithm should also be computed.}
  \item{returnall}{A boolean indicating whether the function returns useful variables or just the denoised datapoints.}
}
\details{
  Essentially, this function applies \code{J} times the modified lifting algorithm (that can be found in \code{\link{fwtnpperm}}), and removes the noise from all sets of detail coefficients by using empirical Bayes shrinkage (of package 'EbayesThresh'). Inverting (by means of the function \code{invtnp} of the package 'adlift') each transform consequently results in \code{J} estimates of the (unknown) true signal. The average of these estimators is our proposed estimator. The functions that appear as choices for \code{Pred} can be found in the package 'adlift'.
}
\value{
  \item{vec}{Matrix whose rows give the trajectories to be used by the nondecimated lifting algorithm.}
  \item{ghatnat}{Vector that gives the estimated true signal given by denoising using the lifting scheme that establishes its own order for removing the points (but with the same specification for prediction stage and neighbourhood as the modified algorithm), rather than a randomly generated order.}
  \item{aveghat}{Estimated signal, obtained as the average of the individual estimates from the random trajectory runs.}
}
\references{See the paper 'A "nondecimated" lifting transform.' by Knight, M.I. and Nason, G.P. (2009) for further details.}
\author{Marina Knight (marina.knight@bristol.ac.uk)}
\note{Use this function together with the "adlift" and "EbayesThresh" packages available from CRAN.}
\seealso{\code{\link{denoiseperm}}, \code{\link{fwtnpperm}}, \code{\link{fwtnpperm}}, and also \code{invtnp} of package 'adlift'}
\examples{
#construct the grid
x<-runif(256) 

#construct the true, normally unknown, signal
g<-make.signal2("blocks",x=x) 

#generate noise with mean 0 and signal-to-noise ratio 5
noise<-rnorm(256,mean=0,sd=sqrt(var(g))/5)

#generate a noisy version of g 
f<-g+noise 

#decide on a number of random trajectories to be used (below J=100, in paper J=20,30), and apply the nondecimated lifting transform to the noisy signal (x,f) 
#below we apply the modified lifting transform J times, each time following a different path, and using adaptive prediction with neighbourhoods of size 2
#in closest configuration; all details are then thresholded using posterior medians and the algorithms inverted
#the aggregate estimator of g proposed by our method is found in out$aveghat
out<-nlt(x,f,J=10,Pred=AdaptPred,neighbours=2,closest=TRUE,intercept=TRUE,nkeep=2,trule="median") 
}
\keyword{nonparametric}
