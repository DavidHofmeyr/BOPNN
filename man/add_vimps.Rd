\name{add_vimps}
\alias{add_vimps}
\title{Compute variable importance scores}
\description{
  Compute variable importance scores for a pre-fitted BOPNN model
}
\usage{
  add_vimps(model)
}
\arguments{
  \item{model}{A BOPNN classification model, output from function BOPNN()}
}
\value{
  The same as input but with entry $vimps updated.
}
\references{
  Hofmeyr, D.P. (2025) "Bags of Projected Nearest Neighbours: Competitors to Random Forests?", ArXiv preprint.
}
\keyword{file}