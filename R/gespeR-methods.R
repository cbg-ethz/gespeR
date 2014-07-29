#' Deconvolute off-target confounded RNAi knockdown phenotypes.
#' 
#' This generic function can handle different types of inputs for phenotype data and target relation matrices. It either reads from two .rds files,
#' uses a data.frame and matrix, respectively, or re-fits the model of a \linkS4class{gespeR} object. The actual model is a regularised linear
#' regression model using an elastic net penalty: \eqn{X = X\beta + \epsilon}. The model can be tuned using the regularisation parameters .. (text see deepSNV!)
#'
#' @rdname gespeR-methods
#' @author Fabian Schmich
#' @exportMethod gespeR
#'
#' @param phenotypes The siRNA-spefic phenotypes
#' @param target.relations The siRNA-to-gene target relations
#' @param mode The mode of covariate selectino ("cv" or "stability")
#' @param ... Additional arguments
#' @return A \linkS4class{gespeR} object
setGeneric("gespeR",
           function(phenotypes, target.relations, ...) {
             standardGeneric("gespeR")
           })

#' @rdname gespeR-methods
setMethod("gespeR",
          signature = signature(phenotypes="Phenotypes", target.relations="TargetRelations"),
          function(phenotypes, target.relations, mode=c("cv", "stability"), ...){
            mode = match.arg(mode)
            if (!target.relations@loaded) target.relations <- loadMatrix(target.relations)
            model <- switch(mode,
                            "cv"=.gespeR.cv(targets=target.relations@values, SSP=scores(phenotypes), alpha=alpha, ...),
                            "stability"=.gespeR.stability(targets=target.relations@values, SSP=scores(phenotypes), ...)
            )            
            gsp <- switch(mode,
                          cv=model$coefficients[-1][which(model$coefficients[-1] != 0)],
                          stability=model$coefficients
            )
            names(gsp) <- switch(mode,
                                 cv=names(gsp),
                                 stability=names(model$stability$selection)
            )
            target.relations <- unloadMatrix(target.relations)
            new("gespeR",
                SSP=phenotypes,
                GSP=Phenotype(phenotype=gsp, ids=names(gsp), type="GSP"),
                target.relations=target.relations,
                model=model,
                is.fitted=TRUE)
            
          }
)