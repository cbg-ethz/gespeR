#' Class definition for Phenotypes
#' 
#' Class used to represent various types of phenotypes, e.g. from 
#' siRNA-specific knockdowns (SSP) or estimated gene-specific phenotypes (GSP).
#' 
#' @author Fabian Schmich
#' @name Phenotypes
#' @rdname Phenotypes-class
#' @seealso \code{\link{Phenotypes}}
#' 
#' @example inst/example/gespeR-example.R
#' @exportClass Phenotypes
#' 
#' @slot type The indicator of represented phenotypes (i.e., "SSP" or "GSP")
#' @slot ids The phenotype identifiers (i.e., siRNA or gene ids)
#' @slot values The phenotypic values
setClass(Class="Phenotypes",
         representation=representation(
           type="character",
           ids="character",
           values="numeric"
         ),
         validity=function(object) {
           if(length(object@ids) != length(object@values)) {
             stop("Number of values does not correspond to number of IDs.")
           }
           return(TRUE)
         }
)
