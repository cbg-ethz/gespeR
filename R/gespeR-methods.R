#' Deconvolute off-target confounded RNAi knockdown phenotypes.
#' 
#' This generic function can handle different types of inputs for phenotype data and target relation matrices. It either reads from two .rds files,
#' uses a data.frame and matrix, respectively, or re-fits the model of a \code{\linkS4class{gespeR}} object. The actual model is a regularised linear
#' regression model using an elastic net penalty: \eqn{X = X\beta + \epsilon}. The model can be tuned using the regularisation parameters .. (text see deepSNV!)
#' 
#' @author Fabian Schmich
#' @rdname gespeR-methods
#' @aliases gespeR
#' @exportMethod gespeR
#'
#' @param phenotypes The siRNA-spefic phenotypes
#' @param target.relations The siRNA-to-gene target relations
#' @param mode The mode of covariate selectino ("cv" or "stability")
#' @param ... Additional arguments
#' @return A \code{\linkS4class{gespeR}} object
setGeneric("gespeR",
           function(phenotypes, target.relations, ...) {
             standardGeneric("gespeR")
           })

#' @rdname gespeR-methods
setMethod("gespeR",
          signature = signature(phenotypes="Phenotypes", target.relations="TargetRelations"),
          function(phenotypes, target.relations, mode=c("cv", "stability"), ...){
            mode = match.arg(mode)
            if (!target.relations@is.loaded) target.relations <- loadValues(target.relations)
            xy <- join(targets=target.relations, phenotypes=phenotypes)
            model <- switch(mode,
                            "cv"=.gespeR.cv(targets=values(xy$targets), SSP=scores(xy$phenotypes), ...),
                            "stability"=.gespeR.stability(targets=values(xy$targets), SSP=scores(xy$phenotypes), ...)
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

#' Constructor for Phenotypes objects
#' 
#' This generic function can handle different types of inputs for phenotype data. It either takes a numeric vector or reads the values from a .txt file.
#' 
#' @author Fabian Schmich
#' @rdname Phenotypes-methods
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
#' @return A \code{\link{Phenotypes}} object 
setGeneric(name="Phenotypes", 
           def=function(phenotypes, ...) { 
             standardGeneric("Phenotypes")
           })

#' @rdname Phenotypes-methods
#' @param sep The separator string
#' @param col.id Column number for the identifier
#' @param col.score Column number for the phenotype score
#' @aliases Phenotypes
setMethod(f="Phenotypes",
          signature=signature(phenotypes="character"),
          function(phenotypes, type=c("SSP", "GSP"), sep="\t", col.id=1, col.score=2) {
            type <- match.arg(type)
            if (!file.exists(phenotypes)) {
              stop(sprintf("File not found: %s", phenotype))
            } else {
              p <- read.delim(phenotypes, sep=sep, stringsAsFactors=F)
            }
            new("Phenotypes", type=type, ids=p[,col.id], values=p[,col.score])
          }
)

#' @rdname Phenotypes-methods
#' @aliases Phenotypes
#' @param channel The cellHTS channel identifier
#' @param sample The cellHTS sample index
setMethod(f="Phenotypes",
          signature=signature(phenotypes="cellHTS"),
          function(phenotypes, channel, sample) {
            new("Phenotypes", type="SSP", ids=featureNames(phenotypes), values=Data(phenotypes)[,,channel])
          }
)


#' @rdname Phenotypes-methods
#' @aliases Phenotypes
setMethod(f="Phenotypes",
          signature=signature(phenotype="numeric"),
          function(phenotypes, ids=NULL, type=c("SSP", "GSP")) {
            type <- match.arg(type)
            if (is.null(ids)) {
              ids <- paste("id", 1:length(phenotypes), sep="_")
            }
            new("Phenotypes", type=type, ids=ids, values=phenotypes)
          }
)



#' scores
#' 
#' Return a named vector of phenotype scores
#' 
#' @author Fabian Schmich
#' @name scores
#' @rdname scores-method
#' @seealso \code{\linkS4class{gespeR}}, \code{\linkS4class{Phenotypes}}
#'  
#' @include gespeR-class.R
#' @exportMethod scores
#' 
#' @param object A \code{\linkS4class{gespeR}} or \code{\linkS4class{Phenotypes}} object
#' @param type The type of phenotype scores (GSP, SSP)
#' @param ... Additional arguments
#' @return A named vector of scores for each phenotype identifier
if (!isGeneric("scores")) {
  setGeneric(name="scores",
           def=function(object, ...) {              
             standardGeneric("scores")
           },
           package="gespeR")
}

#' @rdname scores-method
setMethod(f="scores",
          signature=signature(object="Phenotypes"),
          definition=function(object) {
            r <- object@values
            names(r) <- object@ids
            return(r)
          }
)

#' @rdname scores-method
setMethod(f="scores",
          signature=signature(object="gespeR"),
          definition=function(object, type=c("GSP", "SSP")) {
            type <- match.arg(type)
            phen <- switch(type,
                           "GSP"=object@GSP,
                           "SSP"=object@SSP
            )
            return(scores(phen))
          }
)


#' join
#' 
#' Join a TargetRelations object and a Phenotype object
#' 
#' @author Fabian Schmich
#' @rdname join-methods
#' @export
#' 
#' @param targets A \code{\linkS4class{TargetRelations}} object.
#' @param phenotypes A \code{\linkS4class{Phenotypes}} object.
#' @return List containing the matched targets and phenotypes
setGeneric(name="join", 
           def=function(targets, phenotypes) {
             standardGeneric("join")
           })

#' @author Fabian Schmich
#' @rdname join-methods
setMethod(f="join",
          signature=c(targets="TargetRelations", phenotypes="Phenotypes"),
          def=function(targets, phenotypes) {
            p.siRNAs <- phenotypes@ids
            t.siRNAs <- targets@siRNAs
            isect <- intersect(p.siRNAs, t.siRNAs)
            isect <- sort(isect)
            if (length(isect) != length(p.siRNAs) | length(isect) != length(t.siRNAs)) {
              warning("Phenotype siRNAs and TargetRelations siRNAs do not fully overlap. Stripping to intersection...")
            }
            if (length(isect) == 0) {
              stop("Phenotype siRNAs and TargetRelations siRNAs do not overlap.")
            }
            phenotypes <- phenotypes[match(isect, p.siRNAs)]
            targets <- targets[match(isect, t.siRNAs),]
            return(list(phenotypes=phenotypes, targets=targets))
          }
)


#' annotate.gsp
#' 
#' Query Biomart HGNC symbols for the entrez identifiers of estimated GSPs.
#' Currently, only implemented for species "hsapiens".
#' 
#' @author Fabian Schmich
#' @rdname annotate.gsp-methods
#' @name annotate.gsp
#' @importFrom biomaRt useMart useDataset getBM
#' @export
#' 
#' @param object A \code{\linkS4class{gespeR}} or \code{\linkS4class{Phenotypes}} object
#' @param organism String indicating the biomaRt organism, e.g. "hsapiens"
#' @param ... Additional arguments
#' @return data.frame containing gene identifier, gene symbol and phenotypic score
if (!isGeneric("annotate.gsp")) {
  setGeneric(name="annotate.gsp",
             def=function(object, ...) {              
               standardGeneric("annotate.gsp")
             },
             package="gespeR")
}

#' @rdname annotate.gsp-methods
setMethod(f="annotate.gsp",
          signature=c(object="Phenotypes"),
          def=function(object, organism="hsapiens") {
            if (length(organism) > 1) organism <- organism[1]
            if (object@type != "GSP") stop("Annotation only for GSP phenotypes available.")
            # Set up biomaRt
            mart <- useMart("ensembl")
            ensembl <- useDataset(dataset=paste(organism, "gene", "ensembl", sep="_"), mart=mart, verbose=FALSE)
            # Retrieve results
            result <- data.frame(GeneID=object@ids, Score=object@values)
            result$GeneID <- gsub("ENTREZ", replacement="", result$GeneID, ignore.case=TRUE)
            symbols <- getBM(attributes=c("entrezgene", "hgnc_symbol"), filters="entrezgene", values=result$GeneID, mart=ensembl)
            if (nrow(symbols) == 0) stop("No annotation found for IDs.")
            # Prepare output
            result <- merge(result, symbols, by.x="GeneID", by.y="entrezgene", all.x=TRUE)
            return(data.frame(GeneID=result$GeneID, GeneSymbol=result$hgnc_symbol, Score=result$Score))
          }
)

#' @rdname annotate.gsp-methods
setMethod(f="annotate.gsp",
          signature=c(object="gespeR"),
          def=function(object, organism="hsapiens") {
            if (!object@is.fitted) stop("Fit model first.")
            return(annotate.gsp(object@GSP, organism))
          }
)



#' Constructor for TargetRelations objects
#' 
#' This generic function can handle different types of inputs for TargetRelations data. It either takes ...
#' 
#' @author Fabian Schmich
#' @rdname TargetRelations-methods
#' @export 
#' 
#' @param targets Path to a .rds target relations matrix file or matrix object
#' @return A \code{\linkS4class{TargetRelations}} object 
setGeneric(name="TargetRelations", 
           def=function(targets) { 
             standardGeneric("TargetRelations")
           })

#' @rdname TargetRelations-methods
#' @aliases TargetRelations
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


#' @rdname TargetRelations-methods
#' @aliases TargetRelations
setMethod(f="TargetRelations",
          signature=signature(targets="Matrix"),
          function(targets) {
            new("TargetRelations",
                path=NULL,
                siRNAs=rownames(targets),
                values=targets,
                genes=colnames(targets),
                is.loaded=TRUE
            )  
          }
)

#' Load values of a TargetRelations object
#' 
#' @author Fabian Schmich
#' @rdname loadValues-methods
#' 
#' @seealso \code{\link{unloadValues}}
#' @export loadValues
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @return A \code{\linkS4class{TargetRelations}} object
setGeneric(name="loadValues", 
           def=function(object) {
             standardGeneric("loadValues")
           })

#' @rdname loadValues-methods
setMethod(f="loadValues",
          signature=signature("TargetRelations"),
          function(object) {
            if (!object@is.loaded) {
              if (file.exists(object@path)) {
                new("TargetRelations",
                    is.loaded=TRUE,
                    siRNAs=object@siRNAs,
                    genes=object@genes,
                    path=object@path,
                    values=readRDS(object@path)
                )
              } else {
                stop(sprintf("File not found: %s", object@path))
              }
            } else {
              return(object)
            }
          }
)

#' @rdname loadValues-methods
setMethod(f="loadValues",
          signature=signature("gespeR"),
          function(object) {
            object@target.relations <- loadValues(object@target.relations)
            return(object)
          }
)

#' Unload values of a TargetRelations object
#' 
#' @author Fabian Schmich
#' @rdname unloadValues-methods
#' 
#' @seealso \code{\link{loadValues}}
#' @export unloadValues
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object or \code{\linkS4class{gespeR}} object
#' @param writeValues Indicator, whether to write values
#' @param overwrite Indicator, wheter to overwrite values if file exists at path
#' @param path The path to write out values
#' @param ... Additional arguments
#' @return A \code{\linkS4class{TargetRelations}} object or \code{\linkS4class{gespeR}} object
setGeneric(name="unloadValues", 
           def=function(object, ...) {
             standardGeneric("unloadValues")
           })

#' @rdname unloadValues-methods
setMethod(f="unloadValues",
          signature=signature(object="TargetRelations"),
          function(object, writeValues=TRUE, overwrite=FALSE) {
            if (writeValues) {
              mat.written <- FALSE
              if (!object@is.loaded)
                object <- loadValues(object)
              mat.written <- writeValues(object=object, overwrite=overwrite)
              if (mat.written) {
                new("TargetRelations",
                    is.loaded=FALSE,
                    siRNAs=object@siRNAs,
                    genes=object@genes,
                    path=object@path,
                    values=Matrix(0,0,0)
                )
              } else {
                #object <- unloadValues(object, writeValues=FALSE, overwrite=FALSE)
                return(object)
              }               
            } else {
              if (object@is.loaded) {
                new("TargetRelations",
                    is.loaded=FALSE,
                    siRNAs=object@siRNAs,
                    genes=object@genes,
                    path=object@path,
                    values=Matrix(0,0,0)
                )
              } else {
                #object <- unloadValues(object, writeValues=FALSE, overwrite=FALSE)
                return(object)
              }
            }
          }
)

#' @rdname unloadValues-methods
setMethod(f="unloadValues",
          signature=signature(object="gespeR"),
          function(object, writeValues=TRUE, overwrite=FALSE, path=NULL) {
            x <- object@target.relations
            if (!is.null(path)) path(x) <- path
            x <- unloadValues(x, writeValues=writeValues, overwrite=overwrite)
            object@target.relations <- x
            return(object)
          }
)

#' Write down the values of a \code{\link{TargetRelations}} object
#' 
#' @author Fabian Schmich
#' @rdname writeValues-methods
#' 
#' @seealso \code{\link{loadValues}} \code{\link{unloadValues}}
#' @export
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @param overwrite Indicator, whether to overwrite values if file exists at path
#' @param ... Additional arguments
#' @return A \code{\linkS4class{TargetRelations}} object
setGeneric(name="writeValues", 
           def=function(object, ...) {
             standardGeneric("writeValues")
           })

#' @rdname writeValues-methods
setMethod(f="writeValues",
          signature=signature("TargetRelations"),
          function(object, overwrite=FALSE) {
            if (!object@is.loaded)
              object <- loadValues(object)
            if (!is.null(object@path)) {
              if (!overwrite & file.exists(object@path)) {
                warning("Values file already exists on HDD. Set overwrite=TRUE.")                
              } else {
                saveRDS(object@values, file=object@path)
                return(TRUE)
              }
            } else {
              warning("No path defined. Use path()")
            }
            return(FALSE)
          }
)


#' Set the path of a A \code{\linkS4class{TargetRelations}} object object
#' 
#' @author Fabian Schmich
#' @export 
#' @rdname path-methods
#' 
#' @param object A \code{\link{TargetRelations}} object
#' @param value A string defining the path
#' @return A \code{\linkS4class{TargetRelations}} object with set path
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


#' Get GSPs
#'
#' @author Fabian Schmich
#' @rdname gsp-methods
#' @export
#' @param object A \code{\linkS4class{gespeR}} object 
#' @return A \code{\linkS4class{Phenotypes}} object of GSPs
setGeneric(name="gsp", def=function(object) standardGeneric("gsp"))
#' @rdname gsp-methods
setMethod(f="gsp",
          signature=signature(object="gespeR"),
          function(object) object@GSP
         )

#' Get SSPs
#'
#' @author Fabian Schmich
#' @rdname ssp-methods
#' @export
#' @param object A \code{\linkS4class{gespeR}} object 
#' @return A \code{\linkS4class{Phenotypes}} object of SSPs
setGeneric(name="ssp", def=function(object) standardGeneric("ssp"))
#' @rdname ssp-methods
setMethod(f="ssp",
          signature=signature(object="gespeR"),
          function(object) object@SSP
)

#' Get stability
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

#' Get values
#'
#' @author Fabian Schmich
#' @rdname values-methods
#' @export
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @return A \code{\linkS4class{Matrix}} object of TargetRelations values
setGeneric(name="values", def=function(object) standardGeneric("values"))
#' @rdname values-methods
setMethod(f="values",
          signature=signature(object="TargetRelations"),
          function(object) {
            if (!object@is.loaded) object <- loadValues(object)
            return(object@values)
          }
)

#' Get target.relations
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

