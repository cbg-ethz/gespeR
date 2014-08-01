#' Package: Gene-Specific Phenotype EstimatoR
#' 
#' This package provides a model to deconvolute off-target confounded RNAi knockdown phentypes, and
#' methods to investigate concordance between ranked lists of (estimated) phenotypes. The regularized
#' linear regression model can be fitted using two different strategies. (a) Cross-validation over
#' regularization parameters optimising the mean-squared-error of the model and (b) stability selection
#' of covariates (genes) based on a method by Nicolai Meinshausen et al.
#' 
#' @name gespeR-package
#' @author Fabian Schmich | Computational Biology Group, ETH ZURICH | \email{fabian.schmich@@bsse.ethz.ch}
#' @seealso \code{\link{gespeR}}
#' @example inst/example/gespeR-example.R
#' @docType package
#' @keywords package
#' @references Fabian Schmich et. al,  Deconvoluting Off-Target Confounded RNA Interference Screens (2014).
#' @import methods Matrix glmnet
#' @importFrom graphics plot

NA

#TODO: add data sets