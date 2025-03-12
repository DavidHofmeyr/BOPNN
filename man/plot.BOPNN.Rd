\name{plot.BOPNN}
\alias{plot.BOPNN}
\alias{plot}
\title{Plot method for class BOPNN}
\description{
  Plots data projected into ensemble discriminant subspace from a pre-fitted BOPNN model
}
\usage{
  \method{plot}{BOPNN}(model, new_data = NULL, labels = NULL, dims = c(1, 2), ...)
}
\arguments{
  \item{model}{BOPNN model, output from function BOPNN()}
  \item{new_data}{(optional) data frame with the same variables used to fit model. If omitted then the training data are shown.}
  \item{labels}{(optional) class labels, potentially predictions output from predict(model, new_data). If omitted all points shown in black, otherwise differentiated by labels.}
  \item{dims}{Integer vector of length two, the dimensions of the discriminant subspace to be shown. Default is the first two.}
  \item{...}{(optional) further arguments passed to or from other methods.}
}
\examples{
  ### See example in help(BOPNN-package)
}
\keyword{file}