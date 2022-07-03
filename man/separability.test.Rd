\name{separability.test}
\alias{separability.test}
\title{
Separability test for spatio-temporal point patterns
}
\description{
Performs a separability test of the first-order intensity function based on a Fisher Monte Carlo test of cell counts.
}
\usage{
separability.test(X, t,
                  nx = NULL, ny = NULL, nt = NULL,
                  nperm = 1000)
}
\arguments{
  \item{X}{A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations}
  \item{t}{A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.}
 \item{nx, ny, nt}{Numbers of quadrats in the \eqn{x,y} and \eqn{t} directions.}
 \item{nperm}{An integer specifying the number of replicates used in the Monte Carlo test.}
}
\details{
This function performs a basic test of the separability hypothesis in a manner similar to independence test in two-way contingency tables. The test is conditional on the observed number of points. It considers a regular division of the interval \eqn{T} into disjoint sub-intervals \eqn{T_1,...,T_{n_t}} and similarly a division of the window \eqn{W} into disjoint subsets \eqn{W1, ..., W{n_x \timex n_y}}. Then the function computes Fisher's test statistic and get a p-value based on Monte Carlos approximation.
}
\value{
A list with class "htest" containing the following components:
\item{p.value}{the approximate p-value of the test.}
\item{method}{the character string "Separability test based on Fisher's for counting data".}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Ghorbani et al. (2021) Testing the first-order separability hypothesis for spatio-temporal point patterns, \emph{Computational Statistics & Data Analysis,} Vol. 161, p. 107245}
}
\author{
Jonatan A. González
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
X <- rpoispp(1400)
t <- rbeta(X$n, 1,4,0.8)
stIntensity <- dens.st.bin(X,t, at = "bins")
plot(as.imlist(stIntensity[13:16]), main = 'Non-separable Example')
}
}
