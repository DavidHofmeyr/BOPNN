\name{BOPNN}
\alias{BOPNN}
\title{Bags of Projected Nearest Neighbours}
\description{
  Fits a classification model based on bootstrap aggregation of nearest neighbours models fit to discriminant subspaces
}
\usage{
  BOPNN(data_df, nn = 3, ndim = NULL, dim_prop1 = .5, ndim_samp = NULL, dim_prop0 = min(.75, 5/sqrt(ncol(data_df))), n_prop = 0.63, nbag = 100, coarseness = 6, vimp = FALSE, eps = 0, cores = 1, OOB = FALSE, case_weights = NULL, plotMat = FALSE)
}
\arguments{
  \item{data_df}{Data frame containing at least two covariates and class labels. The entry in data_df for the labels MUST be named "y".}
  \item{nn}{Integer number of nearest neighbours used in each model. Default is 3.}
  \item{ndim}{(optional) Integer number of dimensions for each subspace. Preferably specify the proportion of randomly sampled variables using dim_prop1.}
  \item{dim_prop1}{Proportion (i.e. lying in (0, 1)) of random subspace used as discriminant subspace. Default is 0.5}
  \item{ndim_samp}{(optional) Integer number of variables to randomly sample for each model. Preferable specify the proportion of total covariates using dim_prop0.}
  \item{dim_prop0}{Proportion (i.e. lying in (0, 1)) of total covariates sampled for each model. Default is minimum of 0.75 and 5/sqrt(total covariates)}
  \item{n_prop}{Proportion of sample used in each bootstrap sample. Default is 0.63 sampled without replacement. To (simply) sample with replacement set case_weights all equal to one.}
  \item{nbag}{Integer number of models in ensemble. Default is 100}
  \item{coarseness}{Integer. All factor and integer variables are automatically one-hot-encoded, as are those with at most coarseness unique values. Default is 6.}
  \item{vimp}{Logical. Should variable importance scores be computed. Default is FALSE. Note that these can be added later using the function add_vimps()}
  \item{eps}{Positive numeric. Error tolerance for approximate nearest neighbour search. Default is 0 (exact search).}
  \item{cores}{Integer number of cores to use for parallel computing. Default is 1 (no parallelisation).}
  \item{OOB}{Logical. Should Out-Of-Bag estimate for performance be computed.}
  \item{case_weights}{Non-negative numeric vector of length equal to number of data. Samples with higher case_weights are sampled with appropriately higher probability. Note that if case_weights are supplied then sampling is with replacement (even if they are all set equal), otherwise sampling is without replacement.}
  \item{plotMat}{Logical. Should an ensemble discriminant projection be computed for plotting purposes. Note that this can also be added later using the function add_plotMat().}
}
\value{
  A named list of class BOPNN with fields (at most, depending on settings)
  \item{$mods}{A list containing all the individual models in the ensemble.}
  \item{$nn}{The number of nearest neighbours used. Used for prediction}
  \item{$convert}{The indices of covariates one-hot-encoded. Used for prediction}
  \item{$idxs}{Mapping of covariates to columns of one-hot-encoded data. Used for prediction}
  \item{$ylevs}{Names of classes}
  \item{$y}{Integer encoded class labels}
  \item{$X}{One-hot-encoded covariates}
  \item{$x0}{Original covariates}
  \item{$scales}{Standardisation coefficients for one-hot-encoded data. Used for prediction}
  \item{$vimp}{Variable importance scores.}
  \item{$dr_all}{Ensemble discriminant projection. If plotMat set to TRUE}
  \item{$probs_oob}{Out-of-bag class probabilities. If OOB set to TRUE}
  \item{$predicted.oob}{Out-of-bag class predictions. If OOB set to TRUE}
  \item{$oob.uncertainty}{Entropy based estimates for out-of-bag prediction uncertainty. If OOB set to TRUE}
  \item{$ooberr}{Out-of-bag estimate of classification error. If OOB is set ot TRUE}\item{$probs_oob}{Out-of-bag class probabilities. If OOB set to TRUE}
  \item{$predicted.oob}{Out-of-bag class predictions. If OOB set to TRUE}
  \item{$oob.uncertainty}{Entropy based estimates for out-of-bag prediction uncertainty. If OOB set to TRUE}
  \item{$ooberr}{Out-of-bag estimate of classification error. If OOB is set ot TRUE}
}
\references{
  Hofmeyr, D.P. (2025) "Bags of Projected Nearest Neighbours: Competitors to Random Forests?", ArXiv preprint.
}
\examples{
  ### See example at help(BOPNN-package)
}
\keyword{file}