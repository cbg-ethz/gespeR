#' Package: Gene-Specific Phenotype EstimatoR
#' 
#' This package provides a model to deconvolute off-target confounded RNAi knockdown phentypes, and
#' methods to investigate concordance between ranked lists of (estimated) phenotypes. The regularized
#' linear regression model can be fitted using two different strategies. (a) Cross-validation over
#' regularization parameters optimising the mean-squared-error of the model and (b) stability selection
#' of covariates (genes) based on a method by Nicolai Meinshausen et al.
#' @name gespeR-package
#' @rdname gespeR-package
#' @aliases gespeRpkg
#' @author Fabian Schmich | Computational Biology Group, ETH ZURICH | \email{fabian.schmich@@bsse.ethz.ch}
#' @example inst/example/gespeR-example.R
#' @docType package
#' @keywords package
#' @references Fabian Schmich et. al, Deconvoluting Off-Target Confounded RNA Interference Screens (2014).
#' @import methods Matrix glmnet
#' @importFrom graphics plot
#' @seealso \code{\link{gespeR}}
NA

#' Example phenotype and target relations data A, B, C, and D
#' 
#' The data set contains simulated data for four screens. Each screen consists of a phenotype vector
#' and target relations between siRNAs and genes, i.e. which siRNA binds which genes (on- and off-targets).
#' @docType data
#' @examples 
#' pheno.a <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package="gespeR"),
#' type = "SSP", col.id = 1, col.score = 2)
#' targets.a <- TargetRelations(system.file("extdata", "TR_screen_A.rds", package="gespeR"))
#' @name simData
NA

#' Example fits for phenotypes from example screening data A, B, C and D
#' 
#' The data set contains four fitted gespeR models using stability selection
#' @docType data
#' @examples data(stabilityfits)
#' @name stabilityfits
NA