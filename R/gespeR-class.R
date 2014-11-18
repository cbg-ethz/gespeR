#' gespeR
#' 
#' Class that represents a gespeR model. It contains a SSP \code{\linkS4class{Phenotypes}} and \code{\linkS4class{TargetRelations}} 
#' representing a siRNA knockdown experiment. When the model is fitted, it additionaly contains estimated GSP \code{\linkS4class{Phenotypes}}.
#' 
#' @author Fabian Schmich
#' @name gespeR-class
#' @rdname gespeR-class
#' @aliases gespeR
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
#' 
#' @seealso \code{\link{gespeR-package}}
#' @seealso \code{\link{plot.gespeR}}
#' @seealso \code{\link{gsp}}
#' @seealso \code{\link{ssp}}
#' @seealso \code{\link{scores}}
#' @seealso \code{\link{stability}}
#' @seealso \code{\link{target.relations}}
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

#' @rdname gespeR-class
#' @exportMethod gespeR
#'
#' @param phenotypes The siRNA-spefic phenotypes
#' @param target.relations The siRNA-to-gene target relations
#' @param mode The mode of covariate selectino ("cv" or "stability")
#' @param alpha The \code{\link{glmnet}} mixing parameter
#' @param nbootstrap The number of bootstrap samples
#' @param fraction The fraction for each bootstrap sample
#' @param threshold The selection threshold
#' @param EV The expected value of wrongly selected elements
#' @param weakness The weakness parameter for randomised lasso
#' @param ncores The number of cores for parallel computation
#' @param ... Additional arguments
#' @return A \code{\linkS4class{gespeR}} object
setGeneric("gespeR",
           function(phenotypes, target.relations, ...) {
             standardGeneric("gespeR")
           })
#' @rdname gespeR-class
setMethod("gespeR",
          signature = signature(phenotypes="Phenotypes", target.relations="TargetRelations"),
          function(phenotypes, target.relations, 
                   mode=c("cv", "stability"), 
                   alpha = 0.5, 
                   nbootstrap = 100, 
                   fraction = 0.67, 
                   threshold = 0.9, 
                   EV = 1, 
                   weakness = 0.8, 
                   ncores = 1,
                   ...){
            mode = match.arg(mode)
            if (!target.relations@is.loaded) target.relations <- loadValues(target.relations)
            phenotypes <- na.rem(phenotypes)
            xy <- join(targets=target.relations, phenotypes=phenotypes)
            model <- switch(mode,
                            "cv"=.gespeR.cv(targets = values(xy$targets), SSP = scores(xy$phenotypes),
                                            alpha = alpha,
                                            ncores = ncores,
                                            ...),
                            "stability"=.gespeR.stability(targets = values(xy$targets), SSP = scores(xy$phenotypes), 
                                                          nbootstrap = nbootstrap,
                                                          fraction = fraction,
                                                          threshold = threshold,
                                                          EV = EV,
                                                          weakness = weakness,
                                                          ncores = ncores,
                                                          ...)
            )            
            gsp <- switch(mode,
                          cv=model$coefficients[-1][which(model$coefficients[-1] != 0)],
                          stability=model$coefficients
            )
            names(gsp) <- switch(mode,
                                 cv=names(gsp),
                                 stability=names(model$stability$selection)
            )
            #             target.relations <- unloadValues(target.relations)
            new("gespeR",
                SSP=xy$phenotypes,
                GSP=Phenotypes(phenotype=gsp, ids=names(gsp), type="GSP"),
                target.relations=xy$targets,
                model=model,
                is.fitted=TRUE)
            
          }
)
#' @rdname gespeR-class
setMethod("gespeR",
          signature = signature(phenotypes="numeric", target.relations="Matrix"),
          function(phenotypes, target.relations, ...){
            phenotypes <- Phenotypes(phenotypes=phenotypes, ids=rownames(phenotypes), type="SSP")
            target.relations <- TargetRelations(target.relations)
            gespeR(phenotypes, target.relations, ...)
          }
)

#' stability
#' 
#' Retrieve a \code{\linkS4class{Phenotypes}} object with stability values from a
#' \code{\linkS4class{gespeR}} object.
#'
#' @author Fabian Schmich
#' @rdname stability-methods
#' @export
#' @param object A \code{\linkS4class{gespeR}} object  
#' @return A \code{\linkS4class{Phenotypes}} object of SSPs
setGeneric(name="stability", def=function(object) standardGeneric("stability"))
#' @rdname stability-methods
setMethod(f="stability",
          signature=signature(object="gespeR"),
          function(object) {
            if (object@is.fitted) {
              if(object@model$type == "stability") {
                return(object@model$stability)
              } else {
                warning("gespeR model not fitted in mode=stability")
              }
            } else {
              warning("gespeR model not fitted")
            }
          }
)

#' Plot method for \code{\linkS4class{gespeR}} objects
#' 
#' @author Fabian Schmich
#' @export
#' @method plot gespeR
#' 
#' @param x A \code{\linkS4class{gespeR}} object
#' @param ... Additional paramters for plot
plot.gespeR <- function(x, ...) {
  if (x@is.fitted) {
    switch(x@model$type,
           "cv"=hist(scores(x, "GSP"), xlab="Scores", main="Gene-Specific Phenotypes", ...),
           "stability"=plot(x@model$coefficients, x@model$stability$frequency[x@model$stability$selection], xlab="Scores", ylab="Stability", main="Gene-Specific Phenotypes", ...)
    )
  } else {
    hist(scores(x, "SSP"), xlab="Scores", main="siRNA-Specific Phenotypes", ...)
  }
}

#' target.relations
#' 
#' Retrieve siRNA-to-gene target relations from a \code{\linkS4class{gespeR}} object.
#'
#' @author Fabian Schmich
#' @rdname target.relations-methods
#' @export 
#' @param object A \code{\linkS4class{gespeR}} object 
#' @return A \code{\linkS4class{TargetRelations}} object 
setGeneric(name="target.relations", def=function(object) standardGeneric("target.relations"))
#' @rdname target.relations-methods
setMethod(f="target.relations",
          signature=signature(object="gespeR"),
          function(object) {
            return(object@target.relations)
          }
)
