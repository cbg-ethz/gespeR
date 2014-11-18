#' @export
if (!isGeneric("summary")) {
  setGeneric(name="summary",
             def=function(object, ...) {
               standardGeneric("summary")
             },
             package="gespeR"
  )  
}
setMethod(f="summary",
          signature=signature(object="gespeR"),
          definition=function(object) {
            cat("gespeR model\n")
            cat("is fitted:", ifelse(object@is.fitted, "YES", "NO"), "\n")
            cat("targets loaded: ", ifelse(object@target.relations@is.loaded, "YES", "NO"), "\n")
            if (object@is.fitted) {
              cat("type:", object@model$type, "\n")
              if(object@model$type == "stability") {
                cat("EV:", object@model$stability$EV, "\n")
                cat("threshold:", object@model$stability$threshold, "\n")
                cat("nbootstrap:", object@model$stability$nbootstrap, "\n")
              }
              if (object@model$type == "cv") {
                cat("alpha:", object@model$cv$alpha, "\n")
              }
              cat("selected genes:", length(scores(object, type=ifelse(object@is.fitted, "GSP", "SSP"))), "\n")
            }            
          }
)

setMethod(f="show",
          signature=signature(object="gespeR"),
          definition=function(object) {
            if (object@is.fitted) {
              res <- switch(object@model$type,
                            "cv"=data.frame(ID=object@GSP@ids, GSP=object@GSP@values, Stability=NA),
                            "stability"=data.frame(ID=object@GSP@ids, GSP=object@GSP@values, Stability=object@model$stability$frequency[object@model$stability$selection])
              )
              res <- res[order(abs(res$GSP), decreasing=T),]
              rownames(res) <- NULL
              if (nrow(res) > 5) {
                print(head(res), 5)
                cat("...\n")            
              } else {
                print(res)
              }
            } else {
              cat("Model not fitted.\n")
            }
          }
)
setMethod("show",
          signature=signature("Phenotypes"),
          function(object) {
            cat(sprintf("%d %s Phenotypes\n\n", length(object@ids), object@type))
            if (length(object@ids) > 5) {
              print(data.frame(IDs=head(object@ids, 5), Scores=head(object@values, 5), row.names = NULL))
              cat("...\n")            
            } else {
              print(data.frame(IDs=object@ids, Scores=object@values))
            }
          }
)
setMethod("show",
          signature=signature("TargetRelations"),
          function(object) {            
            if(object@is.loaded) {
              cat(sprintf("%d x %d siRNA-to-gene relations.\n", nrow(object@values), ncol(object@values)))
              print(object@values[1:min(nrow(object@values), 10), 1:min(ncol(object@values), 5)])
              cat("...\n")
            } else {
              cat("Values not loaded:", object@path, "\n")
            }
          }
)

#' Subsetting for Phenotype objects.
#' 
#' @author Fabian Schmich
#'
#' @param x A \code{\linkS4class{Phenotypes}} object
#' @param i The subsetting indices
#' @param j Subsetting indices (not used)
#' @param drop Drop Redundant Extent Information
#' @param ... Additional parameters
#' @return A \code{\linkS4class{Phenotypes}} object
setMethod(f="[",
          signature=signature(x="Phenotypes", i="ANY"),
          def=function(x, i, ...){
            #             if (any(is.na(i)))
            #               stop("subscript 'i' contains NA")      
            result = new("Phenotypes", 
                         type=slot(x, "type"),
                         ids=slot(x, "ids")[i, drop=FALSE],
                         values=slot(x, "values")[i, drop=FALSE]
            )
            return(result)
          }
)

#' Subsetting for TargetRelations objects.
#' 
#' @author Fabian Schmich
#'
#' @param x A \code{\linkS4class{TargetRelations}} object
#' @param i The row subsetting indices (siRNAs)
#' @param j The column subsetting indeces (genes)
#' @param drop Drop Redundant Extent Information
#' @param ... Additional parameters
#' @return A \code{\linkS4class{TargetRelations}} object
setMethod(f="[",
          signature=signature(x="TargetRelations", i="ANY", j="ANY"),
          def=function(x, i, j, ...){
            if (missing(i)) i <- 1:nrow(x@values)
            if (missing(j)) j <- 1:ncol(x@values)
            #             if (any(is.na(i)) | any(is.na(j)))
            #               stop("subscript 'i' contains NA")
            is.loaded <- x@is.loaded
            if (!is.loaded) x <- loadValues(x)
            result = new("TargetRelations", 
                         siRNAs=slot(x, "siRNAs")[i, drop=FALSE],
                         genes=slot(x, "genes")[j, drop=FALSE],
                         values=slot(x, "values")[i, j, drop=FALSE], 
                         path=slot(x, "path"),
                         is.loaded=is.loaded
            )
            if (!is.loaded) result <- unloadValues(result)
            return(result)
          }
)