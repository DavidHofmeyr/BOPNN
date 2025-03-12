\name{BOPNN-package}
\alias{BOPNN-package}
\docType{package}
\title{
  Bags of Projected Nearest Neighbours
}
\description{
  Classification using a bootstrap aggregation of nearest neighbours classifiers fit on different discriminant subspaces.
}
\details{
  The main function in the package is BOPNN(), which fits a classification model.
}
\examples{
  ### generate some simple Gaussian classes, with some noise dimensions
  n <- 1000 ## number of observations
  d <- 10 ## number of informative dimensions
  dnoise <- 20 ## number of noise dimensions
  nC <- 5 ## number of classes
  
  ### simulate data
  
  X <- matrix(rnorm(n*d), n, d)
  y <- sample(1:nC, n, replace = TRUE)
  for(c in 1:nC) X[y==c,] <- (X[y==c,] + matrix(rnorm(d), sum(y==c), d, byrow = TRUE))\%*\%(matrix(rnorm(d^2), d, d)*matrix(rexp(d), d, d, byrow = TRUE))
  
  X <- cbind(X, matrix(rnorm(n*dnoise), n, dnoise)*d)
  
  ### fit model
  ### Argument plotMat allows for plotting discriminant projections
  ### Argument vimp allows for computing variable importance
  ### Alternatively these can both be added after fitting using the functions add_plotMat and add_vimps
  model <- BOPNN(data.frame(X, y = y), plotMat = TRUE, vimp = TRUE)
  
  ### plot variable importance (see that it doesn't use the noise dimensions)
  plot(model$vimp)
  
  ### plot first two (ensemble) discriminant projections
  plot(model, labels = y)
  
}
\references{
  Hofmeyr, D.P (2025) "Bags of Projected Nearest Neighbours: Competitors to Random Forests?", ArXiv preprint.
}
\author{
  David P. Hofmeyr
}