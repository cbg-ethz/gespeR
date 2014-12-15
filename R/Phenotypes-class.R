#' Phenotypes
#' 
#' Class used to represent various types of phenotypes, e.g. from 
#' siRNA-specific (SSP) or estimated gene-specific phenotypes (GSP).
#' 
#' @author Fabian Schmich
#' @name Phenotypes-class
#' @rdname Phenotypes-class
#' @aliases Phenotypes
#' 
#' @exportClass Phenotypes
#' 
#' @slot type The type of represented phenotypes (i.e., "SSP" or "GSP")
#' @slot ids The phenotype identifiers (i.e., siRNA or gene ids)
#' @slot values The phenotypic values
#' 
#' @seealso \code{\link{plot.Phenotypes}}
#' @seealso \code{\link{join}}
#' @seealso \code{\link{gsp}}
#' @seealso \code{\link{ssp}}
#' @seealso \code{\link{scores}}
#' @seealso \code{\link{concordance}}
#' 
#' @examples
#' phenos <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
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

#' @rdname Phenotypes-class
#' 
#' @importClassesFrom cellHTS2 cellHTS
#' @importFrom cellHTS2 Data
#' @importFrom Biobase featureNames channelNames
#' 
#' @export Phenotypes
#' 
#' @param phenotypes The phenotypes as numeric vector, path to a .txt file with two columns (1: identifiers, 2: values), or a cellHTS object
#' @param ids The phenotype identifiers
#' @param type The type of phenotype (GSP, SSP)
#' @param ... Additional arguments
#' @return A \code{\linkS4class{Phenotypes}} object 
setGeneric(name = "Phenotypes", 
           def = function(phenotypes, ...) { 
             standardGeneric("Phenotypes")
           })

#' @rdname Phenotypes-class
#' @param sep The separator string
#' @param col.id Column number for the identifier
#' @param col.score Column number for the phenotype score
setMethod(f="Phenotypes",
          signature = signature(phenotypes = "character"),
          function(phenotypes, type = c("SSP", "GSP"), sep = "\t", col.id = 1, col.score = 2) {
            type <- match.arg(type)
            if (!file.exists(phenotypes)) {
              stop(sprintf("File not found: %s", phenotypes))
            } else {
              p <- read.delim(phenotypes, sep=sep, stringsAsFactors = FALSE)
            }
            new("Phenotypes", type = type, ids = p[,col.id], values = p[,col.score])
          }
)

#' @rdname Phenotypes-class
#' @param channel The cellHTS channel identifier
#' @param sample The cellHTS sample index
setMethod(f="Phenotypes",
          signature=signature(phenotypes = "cellHTS"),
          function(phenotypes, channel, sample) {
            new("Phenotypes", type = "SSP", ids = featureNames(phenotypes), values = Data(phenotypes)[,,channel])
          }
)


#' @rdname Phenotypes-class
setMethod(f="Phenotypes",
          signature = signature(phenotype = "numeric"),
          function(phenotypes, ids = NULL, type = c("SSP", "GSP")) {
            type <- match.arg(type)
            if (is.null(ids)) {
              if (!is.null(rownames(phenotypes))) {
                ids <- rownames(phenotypes)
              } else if (!is.null(names(phenotypes))) {
                ids <- names(phenotypes)
              } else {
                ids <- paste("id", 1:length(phenotypes), sep="_")  
              }
            }
            new("Phenotypes", type = type, ids = ids, values = phenotypes)
          }
)

#' Remove NA/Inf values from phenotype vectors
#' 
#' @author Fabian Schmich
#' @rdname na.rem-methods
#' 
#' @export na.rem
#' 
#' @param object A \code{\linkS4class{Phenotypes}} object
#' @return A \code{\linkS4class{Phenotypes}} object without NA scores values
#' @examples
#' phenos <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' na.rem(phenos)
setGeneric(name="na.rem", 
           def=function(object) {
             standardGeneric("na.rem")
           })

#' @rdname na.rem-methods
setMethod(f="na.rem",
          signature=signature("Phenotypes"),
          function(object) {
            object[which(!is.na(scores(object)))]
          }
)

#' Plot method for Phenotype objects
#' 
#' @return NULL
#' @method plot Phenotypes
#' @export
#' @author Fabian Schmich 
#' @param x A \code{\linkS4class{Phenotypes}} object
#' @param ... Additional arguments for plot
#' @return Histogram of scores
plot.Phenotypes <- function(x, ...) {
  hist(x@values, main=sprintf("%s Phenotypes", x@type), xlab="Scores", ...)
}