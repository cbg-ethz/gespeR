#' TargetRelations
#'
#' Class used to represent siRNA-to-gene on- and off-target relations
#' for a knockdown library and a set of genes.
#' 
#' @author Fabian Schmich
#' @name TargetRelations-class
#' @rdname TargetRelations-class
#' @aliases TargetRelations
#' @exportClass TargetRelations
#' 
#' @slot siRNAs The siRNA identifiers
#' @slot genes The gene identifiers (Entrez)
#' @slot path The path to and \code{.rds} \code{\linkS4class{TargetRelations}} file
#' @slot is.loaded An indicator if target relations values are loaded
#' @slot values The quantitative target relation values between siRNAs and genes
#' 
#' @seealso \code{\link{join}}
#' @seealso \code{\link{loadValues}}
#' @seealso \code{\link{unloadValues}}
#' @seealso \code{\link{writeValues}}
#' @seealso \code{\link{values}}
#' @seealso \code{\link{path<-}}
#' 
#' @examples
#' trels <- TargetRelations(readRDS(system.file("extdata", "TR_screen_A.rds", package = "gespeR")))
setClass(Class="TargetRelations",
         representation=representation(
           siRNAs="character",
           genes="character",
           path="character",
           is.loaded="logical",
           values="Matrix"
         ),
         validity=function(object) {
           if (object@is.loaded) {
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

#' @rdname TargetRelations-class
#' @export 
#' 
#' @param targets Path to a .rds target relations matrix file or \code{\link{Matrix}} object
#' @return A \code{\linkS4class{TargetRelations}} object 
setGeneric(name="TargetRelations", 
           def=function(targets) { 
             standardGeneric("TargetRelations")
           })

#' @rdname TargetRelations-class
setMethod(f="TargetRelations",
          signature=signature(targets="character"),
          function(targets) {
            mat <- readRDS(targets)
            new("TargetRelations",
                path=targets,
                siRNAs=rownames(mat),
                values=mat,
                genes=colnames(mat),
                is.loaded=TRUE
            )  
          }
)

#' @rdname TargetRelations-class
setMethod(f="TargetRelations",
          signature=signature(targets="Matrix"),
          function(targets) {
            new("TargetRelations",
                path="",
                siRNAs=rownames(targets),
                values=targets,
                genes=colnames(targets),
                is.loaded=TRUE
            )  
          }
)

#' path
#' 
#' Set the path of a a \code{\linkS4class{TargetRelations}} object object
#' 
#' @author Fabian Schmich
#' @export 
#' @rdname path-methods
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @param value A string defining the path
#' @return A \code{\linkS4class{TargetRelations}} object with set path
#' @examples
#' trels <- TargetRelations(readRDS(system.file("extdata", "TR_screen_A.rds", package = "gespeR")))
#' path(trels) <- "/dev/null"
setGeneric(name="path<-", 
           def=function(object, value) {
             standardGeneric("path<-")
           })

#' @rdname path-methods
setMethod(f="path<-",
          signature=signature(object="TargetRelations", value="character"),
          function(object, value) {
            object@path <- value
            return(object)
          }
)

#' values
#' 
#' Retrieve the numeric siRNA-to-gene target relations from a \code{\linkS4class{TargetRelations}} object
#'
#' @author Fabian Schmich
#' @rdname values-methods
#' @export
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @return A \code{\link{Matrix}} object
setGeneric(name="values", def=function(object) standardGeneric("values"))
#' @rdname values-methods
#' @examples
#' trels <- TargetRelations(readRDS(system.file("extdata", "TR_screen_A.rds", package = "gespeR")))
#' values(trels)[1:5, 1:5]
setMethod(f="values",
          signature=signature(object="TargetRelations"),
          function(object) {
            if (!object@is.loaded) object <- loadValues(object)
            return(object@values)
          }
)