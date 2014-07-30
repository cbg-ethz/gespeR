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
            if (!target.relations@loaded) target.relations <- loadMatrix(target.relations)
            model <- switch(mode,
                            "cv"=.gespeR.cv(targets=target.relations@values, SSP=scores(phenotypes), ...),
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
                GSP=Phenotypes(phenotype=gsp, ids=names(gsp), type="GSP"),
                target.relations=target.relations,
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
#' @exportMethod Phenotypes
#' 
#' @param phenotypes The phenotypes as numeric vector, path to a .txt file with two columns (1: identifiers, 2: values), or a cellHTS object
#' @param ids The phenotype identifiers
#' @param type The type of phenotype (GSP, SSP)
#' @return A \code{\link{Phenotype}} object 
setGeneric(name="Phenotypes", 
           def=function(phenotypes, ...) { 
             standardGeneric("Phenotypes")
           })

#' @rdname Phenotypes-methods
#' @aliases Phenotypes
setMethod(f="Phenotypes",
          signature=signature(phenotypes="character"),
          function(phenotypes, type=c("SSP", "GSP"), sep="\t", col.id=1, col.score=2, ...) {
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
#' @name join
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
              warning("Phenotype siRNAs and TargetRelations siRNAs do not overlap.")
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
#' @name annotate.gsp
#' @rdname annotate.gsp-methods
#' 
#' @importFrom biomaRt useMart useDataset getBM
#' @export
#' 
#' @param object A \code{\linkS4class{gespeR}} or \code{\linkS4class{Phenotypes}} object
#' @param organism String indicating the biomaRt organism, e.g. "hsapiens"
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
#' @exportMethod TargetRelations
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
                loaded=TRUE
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
                loaded=TRUE
            )  
          }
)

#' Load matrix of a TargetRelations object
#' 
#' @author Fabian Schmich
#' @rdname loadMatrix-methods
#' 
#' @seealso \code{\link{unloadMatrix}}
#' @exportMethod loadMatrix
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @return A \code{\linkS4class{TargetRelations}} object
setGeneric(name="loadMatrix", 
           def=function(object) {
             standardGeneric("loadMatrix")
           })

#' @rdname loadMatrix-methods
setMethod(f="loadMatrix",
          signature=signature("TargetRelations"),
          function(object) {
            if (!object@loaded) {
              if (file.exists(object@path)) {
                new("TargetRelations",
                    loaded=TRUE,
                    siRNAs=object@siRNAs,
                    genes=object@genes,
                    path=object@path,
                    values=readRDS(object@path)
                )
              } else {
                stop(sprintf("File not found: %s", object@path))
              }
            }
          }
)

#' Unload matrix of a TargetRelations object
#' 
#' @author Fabian Schmich
#' @rdname unloadMatrix-methods
#' 
#' @seealso \code{\link{loadMatrix}}
#' @exportMethod unloadMatrix
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @return A \code{\linkS4class{TargetRelations}} object
setGeneric(name="unloadMatrix", 
           def=function(object, ...) {
             standardGeneric("unloadMatrix")
           })

#' @rdname unloadMatrix-methods
setMethod(f="unloadMatrix",
          signature=signature(object="TargetRelations"),
          function(object, writeMatrix=TRUE, overwrite=FALSE) {
            if (writeMatrix) {
              mat.written <- FALSE
              if (!object@loaded)
                object <- loadMatrix(object)
              mat.written <- writeMatrix(object=object, overwrite=overwrite)
              if (mat.written) {
                new("TargetRelations",
                    loaded=FALSE,
                    siRNAs=object@siRNAs,
                    genes=object@genes,
                    path=object@path,
                    values=Matrix(0,0,0)
                )
              } else {
                #object <- unloadMatrix(object, writeMatrix=FALSE, overwrite=FALSE)
                return(object)
              }               
            } else {
              if (object@loaded) {
                new("TargetRelations",
                    loaded=FALSE,
                    siRNAs=object@siRNAs,
                    genes=object@genes,
                    path=object@path,
                    values=Matrix(0,0,0)
                )
              } else {
                #object <- unloadMatrix(object, writeMatrix=FALSE, overwrite=FALSE)
                return(object)
              }
            }
          }
)


#' Write down the values of a \code{\link{TargetRelations}} object
#' 
#' @author Fabian Schmich
#' @rdname writeMatrix-methods
#' 
#' @seealso \code{\link{loadMatrix}} \code{\link{unloadMatrix}}
#' @exportMethod writeMatrix
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object
#' @return A \code{\linkS4class{TargetRelations}} object
setGeneric(name="writeMatrix", 
           def=function(object, ...) {
             standardGeneric("writeMatrix")
           })

#' @rdname writeMatrix-methods
setMethod(f="writeMatrix",
          signature=signature("TargetRelations"),
          function(object, overwrite=FALSE) {
            if (!object@loaded)
              object <- loadMatrix(object)
            if (!is.null(object@path)) {
              if (!overwrite & file.exists(object@path)) {
                warning("Matrix file already exists on HDD. Set overwrite=TRUE.")                
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
#' @exportMethod path<-
#' 
#' @param object A \code{\link{TargetRelations}} object
#' @param path A string defining the path
#' @return A \code{\linkS4class{TargetRelations}} object with set path
setGeneric(name="path<-", 
           def=function(object, value) {
             standardGeneric("path<-")
           })

#' @rdname writeMatrix-methods
setMethod(f="path<-",
          signature=signature(object="TargetRelations", value="character"),
          function(object, value) {
            object@path <- value
            return(object)
          }
)




