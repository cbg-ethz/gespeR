#' scores
#' 
#' Return a named vector of phenotype scores
#' 
#' @author Fabian Schmich
#' @name scores
#' @rdname scores-methods
#'  
#' @include gespeR-class.R
#' @exportMethod scores
#' 
#' @param object A \code{\linkS4class{gespeR}} or \code{\linkS4class{Phenotypes}} object
#' @param type The type of phenotype scores (GSP, SSP)
#' @return A named vector of scores for each phenotype identifier
#' 
#' @seealso \code{\linkS4class{gespeR}}
#' @seealso \code{\linkS4class{Phenotypes}}
#' 
#' @examples
#' data(stabilityfits)
#' scores(stabilityfits$A)
if (!isGeneric("scores")) {
  setGeneric(name="scores",
           def=function(object, ...) {              
             standardGeneric("scores")
           },
           package="gespeR")
}

#' @rdname scores-methods
setMethod(f = "scores",
          signature = signature(object = "Phenotypes"),
          definition = function(object) {
            object %>% as.data.frame() %>% tbl_df()
          }
)

#' @rdname scores-methods
setMethod(f="scores",
          signature=signature(object = "gespeR"),
          definition = function(object, type = c("GSP", "SSP")) {
            type <- match.arg(type)
            phen <- switch(type,
                           "GSP" = object@GSP,
                           "SSP" = object@SSP
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
#' @examples
#' phenos <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' trels <- TargetRelations(readRDS(system.file("extdata", "TR_screen_A.rds", package = "gespeR")))
#' phenos <- phenos[1:17]
#' stripped_down <- join(targets = trels, phenotypes = phenos)
setGeneric(name = "join", 
           def = function(targets, phenotypes) {
             standardGeneric("join")
           })

#' @rdname join-methods
setMethod(f = "join",
          signature = c(targets = "TargetRelations", phenotypes = "Phenotypes"),
          def = function(targets, phenotypes) {
            p.siRNAs <- phenotypes@ids
            t.siRNAs <- targets@siRNAs
            isect <- intersect(p.siRNAs, t.siRNAs)
            isect <- sort(isect)
            if (length(isect) != length(p.siRNAs) | length(isect) != length(t.siRNAs)) {
              #message("Phenotype siRNAs and TargetRelations siRNAs do not fully overlap. Stripping to intersection...")
            }
            if (length(isect) == 0) {
              stop("Phenotype siRNAs and TargetRelations siRNAs do not overlap.")
            }
            phenotypes <- phenotypes[match(isect, p.siRNAs)]
            targets <- targets[match(isect, t.siRNAs),]
            return(list(phenotypes = phenotypes, targets = targets))
          }
)


#' Methods for values of \code{\linkS4class{TargetRelations}} objects
#'
#' Load, unload or write to file the values of a \code{\linkS4class{TargetRelations}} object
#' 
#' @author Fabian Schmich
#' @rdname trmatrix-methods
#' 
#' @export
#' 
#' @param object A \code{\linkS4class{TargetRelations}} object or \code{\linkS4class{gespeR}} object
#' @param writeValues Indicator, whether to write values
#' @param overwrite Indicator, wheter to overwrite values if file exists at path
#' @param path The path to write out values
#' @param ... Additional arguments
#' @return A \code{\linkS4class{TargetRelations}} object or \code{\linkS4class{gespeR}} object
#' @examples
#' data(stabilityfits)
#' \dontrun{
#' loadValues(stabilityfits$A)
#' }
setGeneric(name="loadValues", 
           def=function(object) {
             standardGeneric("loadValues")
           })

#' @rdname trmatrix-methods
setMethod(f="loadValues",
          signature=signature("TargetRelations"),
          function(object) {
            if (!object@is.loaded) {
              if (file.exists(object@path)) {
#                 new("TargetRelations",
#                     is.loaded = TRUE,
#                     siRNAs = object@siRNAs,
#                     genes = object@genes,
#                     path = object@path,
#                     values = readRDS(object@path)
#                 )
                TargetRelations(object@path)
              } else {
                stop(sprintf("File not found: %s", object@path))
              }
            } else {
              return(object)
            }
          }
)

#' @rdname trmatrix-methods
setMethod(f="loadValues",
          signature=signature("gespeR"),
          function(object) {
            object@target.relations <- loadValues(object@target.relations)
            return(object)
          }
)


#' @rdname trmatrix-methods
#' @export
setGeneric(name="unloadValues", 
           def=function(object, ...) {
             standardGeneric("unloadValues")
           })

#' @rdname trmatrix-methods
setMethod(f="unloadValues",
          signature=signature(object="TargetRelations"),
          function(object, writeValues = TRUE, overwrite = FALSE, path = NULL) {
            if (!is.null(path)) path(object) <- path
            if (writeValues) {
              mat.written <- FALSE
              if (!object@is.loaded) object <- loadValues(object)
              mat.written <- writeValues(object=object, overwrite=overwrite)
              if (mat.written) {
                new("TargetRelations",
                    is.loaded = FALSE,
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
                    is.loaded = FALSE,
                    siRNAs = object@siRNAs,
                    genes = object@genes,
                    path = object@path,
                    values = Matrix(0,0,0)
                )
              } else {
                warning("Values already unloaded.")
                return(object)
              }
            }
          }
)

#' @rdname trmatrix-methods
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


#' @rdname trmatrix-methods
#' @export
setGeneric(name="writeValues", 
           def=function(object, ...) {
             standardGeneric("writeValues")
           })

#' @rdname trmatrix-methods
setMethod(f="writeValues",
          signature=signature("TargetRelations"),
          function(object, overwrite=FALSE) {
            if (!object@is.loaded)
              object <- loadValues(object)
            if (!is.null(object@path)) {
              if (!overwrite & file.exists(object@path)) {
                warning("Values file already exists on HDD. Set overwrite = TRUE or change path.")                
              } else {
                saveRDS(object@values, file=object@path)
                return(TRUE)
              }
            } else {
              warning("No path defined. Set path for object.")
            }
            return(FALSE)
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
#' @param organism String indicating the biomaRt organism
#' @return data.frame containing gene identifier, gene symbol and phenotypic score
#' @seealso \code{\link{gsp}}
#' @seealso \code{\link{ssp}}
#' @seealso \code{\link{scores}}
#' @examples
#'  data(stabilityfits)
#'  gspA <- gsp(stabilityfits$A)
#' \dontrun{
#'  annotate.gsp(gspA)
#' }
if (!isGeneric("annotate.gsp")) {
  setGeneric(name = "annotate.gsp",
             def = function(object, ...) {              
               standardGeneric("annotate.gsp")
             },
             package = "gespeR")
}

#' @rdname annotate.gsp-methods
setMethod(f = "annotate.gsp",
          signature = c(object="Phenotypes"),
          def = function(object, organism = "hsapiens") {
            if (length(organism) > 1) organism <- organism[1]
            if (object@type != "GSP") stop("Annotation only for GSP phenotypes available.")
            # Set up biomaRt
            mart <- useMart("ensembl")
            ensembl <- useDataset(dataset = paste(organism, "gene", "ensembl", sep = "_"),
                                  mart = mart, 
                                  verbose  =FALSE)
            # Retrieve results
            result <- as.data.frame(object)
            result$GeneID <- gsub("ENTREZ", replacement="", result$GeneID, ignore.case=TRUE)
            symbols <- getBM(attributes = c("entrezgene", "hgnc_symbol"), 
                             filters = "entrezgene", 
                             values = result$ID, 
                             mart = ensembl)
            if (nrow(symbols) == 0) stop("No annotation found for IDs.")
            # Prepare output
            result <- merge(result, symbols, by.x = "ID", by.y = "entrezgene", all.x = TRUE)
            return(data.frame(GeneID = result$GeneID, GeneSymbol = result$hgnc_symbol, Score = result$Score))
          }
)

#' @rdname annotate.gsp-methods
setMethod(f = "annotate.gsp",
          signature=c(object = "gespeR"),
          def=function(object, organism = "hsapiens") {
            if (!object@is.fitted) stop("Fit model first")
            return(annotate.gsp(object@GSP, organism))
          }
)


#' Retrieve GSPs and SSPs from \code{\linkS4class{gespeR}} objects
#'
#' @author Fabian Schmich
#' @rdname gspssp-methods
#' @export
#' @param object A \code{\linkS4class{gespeR}} object 
#' @return A \code{\linkS4class{Phenotypes}} object of GSPs and SSPs, respectively
#' @seealso \code{\link{annotate.gsp}}
#' @seealso \code{\link{scores}}
#' @examples
#' data(stabilityfits)
#' gsp(stabilityfits$A)
#' ssp(stabilityfits$B)
setGeneric(name = "gsp", def = function(object) standardGeneric("gsp"))

#' @rdname gspssp-methods
setMethod(f = "gsp",
          signature = signature(object = "gespeR"),
          function(object) object@GSP
         )


#' @rdname gspssp-methods
#' @inheritParams gsp
#' @export
setGeneric(name = "ssp", def = function(object) standardGeneric("ssp"))

#' @rdname gspssp-methods
setMethod(f = "ssp",
          signature = signature(object = "gespeR"),
          function(object) object@SSP
)



#' values
#' 
#' Retrieve the numeric values from a \code{\linkS4class{TargetRelations}} 
#' or \code{\linkS4class{Phenotypes}} object
#'
#' @author Fabian Schmich
#' @rdname values-methods
#' @export
#' @param object A \code{\linkS4class{TargetRelations}} 
#' or \code{\linkS4class{Phenotypes}} object
#' @return A \code{\link{Matrix}} object
setGeneric(name = "values", def=function(object) standardGeneric("values"))

#' @rdname values-methods
#' @examples
#' trels <- TargetRelations(readRDS(system.file("extdata", "TR_screen_A.rds", package = "gespeR")))
#' values(trels)[1:5, 1:5]
setMethod(f = "values",
          signature = signature(object = "TargetRelations"),
          function(object) {
            if (!object@is.loaded) object <- loadValues(object)
            return(object@values)
          }
)
#' @rdname values-methods
#' @examples
#' phenos <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' values(phenos)
setMethod(f = "values",
          signature = signature(object = "Phenotypes"),
          function(object) {
            return(object@values)
          }
)