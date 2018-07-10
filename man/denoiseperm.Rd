\name{denoiseperm}
\alias{denoiseperm}
%\alias{denoiseperm2}
\title{Denoise a signal using the modified lifting transform and empirical Bayes thresholding}
\description{
  Denoises an input signal contaminated by noise. First the signal is decomposed using the modified lifting scheme (coded in \code{\link{fwtnpperm}}) using a prespecified order, known as path or trajectory, of point removal. Once the signal is decomposed into wavelet coefficients (or details), these are subjected to an empirical Bayes shrinkage procedure in order to remove the noise, the transform is inverted and an estimate of the noisy signal is obtained.
}
\usage{
denoiseperm(x, f, pred=LinearPred, neigh=1, int=TRUE, clo=FALSE, keep=2, 
rule = "median", per = sample(1:length(x),(length(x)-keep),FALSE),returnall=FALSE)
}
\arguments{
  \item{x}{Vector of any length (not necessarily equally spaced) that gives the grid on which the signal is observed.}
  \item{f}{Vector of the same length as \code{x} that gives the signal values corresponding to the \code{x}-locations. }
  \item{pred}{The type of regression to be used in the prediction step of the modified lifting algorithm. Choices are linear, quadratic or cubic (respectively, \code{LinearPred}, \code{QuadPred} or \code{CubicPred}), or two adaptive procedure which automatically choose the degree used in regression, (\code{AdaptPred} or \code{AdaptNeigh}).}
  \item{neigh}{Number of neighbours to be used in order to construct the neighbourhood of each point that has to be removed. If '\code{clo}=FALSE', this gives the number of neighbours on each side of the removed point.}
  \item{int}{Specifies whether (\code{int}=TRUE) or not (\code{int}=FALSE) an intercept is to be used in the regression curve. For \code{pred}=\code{AdaptPred} or \code{AdaptNeigh}, the algorithm automatically makes this choice.}
  \item{clo}{If (\code{clo}=TRUE) or (\code{clo}=FALSE), then at each step the neighbours are in closest, respectively symmetrical configuration.}
  \item{keep}{Number of scaling points we want at the end of the transform. The usual choice is \code{keep}=2.}
  \item{rule}{The type of Bayesian shrinkage technique, with possible choices posterior median (\code{"median"}) or posterior mean (\code{"mean"}).}
  \item{per}{Vector of length (length(\code{x})-\code{keep}) which gives the order of point removal in the lifting algorithm.}
  \item{returnall}{Indicates whether the function returns useful variables or just the denoised datapoints.}
}
\details{
  Once the modified lifting transform is applied, the wavelet coeficients are divided into artificial levels. The details obtained by means of a lifting scheme have different variances, and will therefore be normalized to have the same variance as the noise. Those normalized details falling into the finest artificial level will be used for estimating the standard deviation of the noise that contaminated the signal. Using this estimate, the normalized details can then be shrunk and un-normalized (using package 'EbayesThresh'), and the transform inverted (using the function \code{invtnp} of package 'adlift') to give an estimate of the signal. The choices for \code{pred} can be found in the package 'adlift'.
}
\value{
 If \env{returnall=FALSE}, the estimate of the function after denoising.  If \env{returnall=TRUE}, a list with components:\cr

  \item{fhat}{Estimated signal after removing the noise.}
  \item{w}{This is the matrix associated to the modified lifting transform.}
  \item{indsd}{Vector giving the standard deviations of the detail and scaling coefficients.}
  \item{al}{List giving the split of points between the artificial levels.}
  \item{sd}{Estimated standard deviation of the noise.}
}
\references{See the paper 'A "nondecimated" lifting transform' by Knight, M.I. and Nason, G.P. (2008) for further details.}
\author{Marina Knight (marina.knight@bristol.ac.uk)}
\note{Use this function together with the "adlift" and "EbayesThresh" packages available from CRAN.}
\seealso{\code{\link{fwtnpperm}}, \code{\link{fwtnpperm}}, and also \code{invtnp} of package 'adlift'}
\examples{
# construct a grid
x<-runif(256)

# construct a true, normally unknown, signal 
g<-make.signal2("bumps",x=x) 

# now generate noise (here with mean 0 and signal-to-noise ratio 3)
noise<-rnorm(256,mean=0,sd=sqrt(var(g))/3) 

# obtain a noisy version of the true signal g
f<-g+noise 

# construct the trajectory which will indicate the order of point removal that will be followed by
# the modified lifting algorithm 
# vec below gives the first (length(x)-keep) entries of a random permutation of (1:length(x))
vec<-sample(1:256,254,FALSE)

# denoise the signal (x,f) by applying the modified lifting transform following the removal order 
# in vec and using adaptive prediction
# and neighbourhoods of size 2 in symmetrical configuration 
# the details are then thresholded using posterior medians and the algorithm inverted
# the proposed estimate of g is given by out$fhat$coeff

out<-denoiseperm(x,f,pred=AdaptPred,neigh=1,int=TRUE,clo=FALSE,keep=2,rule="median",per=vec) 
}
\keyword{nonparametric}
