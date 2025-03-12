\name{add_plotMat}
\alias{add_plotMat}
\title{Compute ensemble projection matrix}
\description{
  Compute ensemble discriminant projection matrix from pre-fitted BOPNN model
}
\usage{
  add_plotMat(model)
}
\arguments{
  \item{model}{A BOPNN classification model, output from function BOPNN()}
}
\value{
  The same as input but with entry $dr_all added.
}
\references{
  Hofmeyr, D.P. (2025) "Bags of Projected Nearest Neighbours: Competitors to Random Forests?", ArXiv preprint.
}
\keyword{file}