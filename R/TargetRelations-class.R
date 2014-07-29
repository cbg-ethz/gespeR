#' Class definition for TargetRelations
#'
#' Class used to represent siRNA-to-gene on- and off-target relations
#' for a knockdown library and a set of genes.
#' 
#' @author Fabian Schmich
#' @name TargetRelations
#' @rdname TargetRelations-class
#' @seealso \code{\link{TargetRelations}}
#' 
#' @example inst/example/gespeR-example.R
#' @exportClass TargetRelations
#' 
#' @slot siRNAs The siRNA identifiers
#' @slot genes The gene identifiers (Entrez)
#' @slot path The path to and \code{.rds} \code{\link{TargetRelations}} file
#' @slot loaded An indicator if target relations values are loaded
#' @slot values The quantitative target relation values between siRNAs and genes
setClass(Class="TargetRelations",
         representation=representation(
           siRNAs="character",
           genes="character",
           path="character",
           loaded="logical",
           values="Matrix"
         ),
         validity=function(object) {
           if (object@loaded) {
             if (!all(dim(object@values) == c(length(object@siRNAs), length(object@genes)))) {
               stop("Incorrect target matrix dimensions.")
             }             
             if(!all(object@genes == colnames(object@values))) {
               stop("Gene IDs do not match with targets matrix.")
             } 
             if(!all(object@siRNAs == rownames(object@values))) {
               stop("siRNA IDs do not match with targets matrix.")
             }
           }
           return(TRUE)
         }
)