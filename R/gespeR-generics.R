#' Summary method for \code{\linkS4class{gespeR}} objects
#' 
#' @author Fabian Schmich
#' @aliases summary, gespeR-method
setMethod(f="summary",
          signature=signature(object="gespeR"),
          definition=function(object) {
            cat("gespeR model\n")
            cat("is.fitted:", ifelse(object@is.fitted, "YES", "NO"), "\n")
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


#' Show method for gespeR objects
#' 
#' @author Fabian Schmich
#' @aliases show, gespeR-method
#' 
#' @param object A \code{\linkS4class{gespeR}} object
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

#' Plot method for \code{\linkS4class{gespeR}} objects
#' 
#' @author Fabian Schmich
#' @export
#' @method plot gespeR
#' 
#' @param x A \code{\linkS4class{gespeR}} object
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



#' Show method for Phenotype objects
#' 
#' @author Fabian Schmich
#' @param object A \code{\link{Phenotype}} object.
setMethod("show",
          signature=signature("Phenotypes"),
          function(object) {
            cat(sprintf("%d %s Phenotypes\n\n", length(object@ids), object@type))
            if (length(object@ids) > 5) {
              print(data.frame(IDs=head(object@ids, 5), Scores=head(object@values, 5)))
              cat("...\n")            
            } else {
              print(data.frame(IDs=object@ids, Scores=object@values))
            }
          }
)

#' Show method for TargetRelations objects
#' 
#' @author Fabian Schmich
#' @param object A \code{\link{TargetRelations}} object.
setMethod("show",
          signature=signature("TargetRelations"),
          function(object) {            
            if(object@loaded) {
              cat(sprintf("%d x %d siRNA-to-gene relations.\n", nrow(object@values), ncol(object@values)))
              print(object@values[1:min(nrow(object@values), 10), 1:min(ncol(object@values), 5)])
              cat("...\n")
            } else {
              cat("Matrix not loaded:", object@path, "\n")
            }
          }
)

#' Plot method for Phenotype objects
#' 
#' @param x A \code{\linkS4class{Phenotypes}} object
#' @return NULL
#' @method plot Phenotypes
#' @export
#' @author Fabian Schmich 
plot.Phenotypes <- function(x, ...) {
  hist(x@values, main=sprintf("%s Phenotypes", x@type), xlab="Scores", ...)
}

#' Subsetting for Phenotype objects.
#' 
#' @author Fabian Schmich
#' @name Extract
#' @name `[`
#' @rdname Extract-methods
#' 
#' @export
#' 
#' @param x A \code{\linkS4class{Phenotypes}} object
#' @param i The subsetting indices
#' @return A \code{\linkS4class{Phenotypes}} object
setMethod(f="[",
          signature=signature(x="Phenotypes"),
          def=function(x, i){
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
#' @name Extract
#' @name `[`
#' @rdname Extract-methods
#' 
#' @export
#' 
#' @param x A \code{\linkS4class{TargetRelations}} object
#' @param i The row subsetting indices (siRNAs)
#' @param j The column subsetting indeces (genes)
#' @return A \code{\linkS4class{TargetRelations}} object
setMethod(f="[",
          signature=signature(x="TargetRelations"),
          def=function(x, i, j){
            if (missing(i)) i <- 1:nrow(x@values)
            if (missing(j)) j <- 1:ncol(x@values)
            #             if (any(is.na(i)) | any(is.na(j)))
            #               stop("subscript 'i' contains NA")
            is.loaded <- x@loaded
            if (!is.loaded) x <- loadMatrix(x)
            result = new("TargetRelations", 
                         siRNAs=slot(x, "siRNAs")[i, drop=FALSE],
                         genes=slot(x, "genes")[j, drop=FALSE],
                         values=slot(x, "values")[i, j, drop=FALSE], 
                         path=slot(x, "path"),
                         loaded=is.loaded
            )
            if (!is.loaded) result <- unloadMatrix(result)
            return(result)
          }
)