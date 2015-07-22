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
#' @slot ids The entity identifiers (i.e., siRNA or gene ids)
#' @slot pnames The phenotype names
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
setClass(Class = "Phenotypes",
         representation = representation(
           type = "character",
           ids = "character",
           pnames = "character",
           values = "Matrix"
         ),
         validity = function(object) {
           if(length(object@ids) != nrow(object@values)) {
             stop("Number of values does not correspond to number of siRNA IDs.")
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
#' @param ids The siRNA/gene identifiers
#' @param pnames The phenotype identifiers
#' @param type The type of phenotype (GSP, SSP)
#' @param ... Additional arguments
#' @return A \code{\linkS4class{Phenotypes}} object 
setGeneric(name = "Phenotypes", 
           def = function(phenotypes, ...) { 
             standardGeneric("Phenotypes")
           })

#' @rdname Phenotypes-class
#' @param sep The separator string
#' @param col.id Column number for the siRNA identifiers
#' @param col.score Column number(s) for the phenotype score
setMethod(f = "Phenotypes",
          signature = signature(phenotypes = "character"),
          function(phenotypes, type = c("SSP", "GSP"), sep = "\t", col.id = 1, col.score = 2) {
            type <- match.arg(type)
            if (!file.exists(phenotypes)) {
              stop(sprintf("File not found: %s", phenotypes))
            } else {
              p <- read.delim(phenotypes, sep=sep, stringsAsFactors = FALSE)
            }
            new("Phenotypes", 
                type = type, 
                ids = p[,col.id],
                pnames = colnames(p[col.score]),
                values = Matrix(as.matrix(p[,col.score], nrow = nrow(p))))
          }
)

#' @rdname Phenotypes-class
#' @param channel The cellHTS channel identifier
#' @param sample The cellHTS sample index
setMethod(f = "Phenotypes",
          signature = signature(phenotypes = "cellHTS"),
          function(phenotypes, channel, sample) {
            new("Phenotypes", 
                type = "SSP", 
                ids = featureNames(phenotypes), 
                pnames = NULL,
                values = Matrix(Data(phenotypes)[,,channel]))
          }
)


#' @rdname Phenotypes-class
setMethod(f = "Phenotypes",
          signature = signature(phenotype = "Matrix"),
          function(phenotypes, ids = NULL, pnames = NULL, type = c("SSP", "GSP")) {
            type <- match.arg(type)
            if (is.null(ids)) {
              if (!is.null(rownames(phenotypes))) {
                ids <- rownames(phenotypes)
              } else {
                ids <- paste("id", 1:nrow(phenotypes), sep="_")  
              }
              if (!is.null(colnames(phenotypes))) {
                pnames <- colnames(phenotypes)
              } else {
                pnames <- paste("phen", 1:ncol(phenotypes), sep = "_")
              }
            }
            new("Phenotypes", 
                type = type, 
                ids = ids,
                pnames = pnames,
                values = phenotypes)
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
setGeneric(name = "na.rem", 
           def = function(object) {
             standardGeneric("na.rem")
           })

#' @rdname na.rem-methods
setMethod(f = "na.rem",
          signature = signature("Phenotypes"),
          function(object) {
            object[apply(object@values, 1, function(x) all(!is.na(x)))]
          }
)

#' Concatenate Phenotypes objects
#' 
#' @author Fabian Schmich
#' @export
#' 
#' @param x A \code{\linkS4class{Phenotypes}} object
#' @param recursive recursive
#' @param ... additional \code{\linkS4class{Phenotypes}} objects
#' @return A concatenated \code{\linkS4class{Phenotypes}} object
#' @examples
#' phenos.a <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' phenos.b <- Phenotypes(system.file("extdata", "Phenotypes_screen_B.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' c(phenos.a, phenos.b)
setMethod("c", 
          signature(x = "Phenotypes"), 
          function(x, ...) {
            elements = list(x, ...)
            t <- unique(sapply(elements, function(x) x@type))
            if (length(t) == 1) {
              if (length(unique(sapply(elements, function(x) length(x@pnames)))) == 1)
                samephens <- all(apply(sapply(elements, function(x) x@pnames) %>% rBind() %>% unique(), 1, function(x) length(unique(x) == 1)))
              else
                samephens <- FALSE
              if (samephens) {
                vals <- do.call("rBind", lapply(elements, function(x) x@values))
                ids <- do.call("c", lapply(elements, function(x) x@ids))
                new("Phenotypes", 
                    type = slot(x, "type"),
                    ids = ids,
                    pnames = slot(x, "pnames"),
                    values = vals)
              } else {
                stop("Unequal phenotypes")
              }
            } else {
              stop("Unequal types")
            }
          }
)

#' Convert Phenotypes object to a data.frame
#' 
#' @author Fabian Schmich
#' 
#' @export
#' 
#' @param x A \code{\linkS4class{Phenotypes}} object
#' @return A data.frame
#' @examples
#' phenos <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' as.data.frame(phenos)
setMethod("as.data.frame", 
          signature(x = "Phenotypes"), 
          function(x) {
            ans <- cbind(data.frame(ID = x@ids), as.data.frame(as.matrix(x@values))) 
            if (!is.null(x@pnames)) colnames(ans) <- c("ID", x@pnames)
            return(ans %>% tbl_df())
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
#' phenos <- Phenotypes(system.file("extdata", "Phenotypes_screen_A.txt", package = "gespeR"),
#' type = "SSP",
#' col.id = 1,
#' col.score = 2)
#' plot(phenos)
plot.Phenotypes <- function(x, main = "") {
  as.data.frame(x) %>% 
    melt(id.vars = "ID") %>% 
    ggplot(aes(x = value)) + 
    geom_histogram() + 
    facet_wrap(~variable) + 
    ggtitle(main) +
#     ggtitle(sprintf("%s Phenotypes", x@type)) +
    xlab("Phenotype") + 
    ylab("Frequency") + 
    theme_bw()
}