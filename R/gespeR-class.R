#' Class definition for gespeR
#' 
#' Class that represents a gespeR model. It includes SSP phenotypes and a TargetRelations object
#' from an siRNA knockdown experiment. When fitted, it additionaly stores the regularized linear regression model
#' and estimated GSP phenotypes.
#' 
#' @author Fabian Schmich
#' @name gespeR-class
#' @rdname gespeR-class
#' @seealso \code{\link{gespeR}}
#' 
#' @include TargetRelations-class.R
#' @include Phenotypes-class.R
#' @example inst/example/gespeR-example.R
#' @exportClass gespeR
#'
#' @slot SSP The observed siRNA-specific phenotypes
#' @slot GSP The deconvoluted gene-specific phenotypes
#' @slot target.relations The siRNA-to-gene target relations, e.g. predicted by TargetScan
#' @slot is.fitted An indicator wheter the gespeR model was fitted
#' @slot model The fitted regularized linear regression model
setClass(Class="gespeR",
         representation=representation(
           SSP="Phenotypes",
           GSP="Phenotypes",
           target.relations="TargetRelations",
           is.fitted="logical",
           model="list"
         ),
         validity=function(object) {
           return(TRUE)
           # TODO: implement validity check
         }
)