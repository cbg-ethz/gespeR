#' Evaluate the concordance between Phenotype objects
#' 
#' Measures include the correlation (rho) between pairs of phenotypes for the same gene, the rank biased overlap (\code{\link{rbo}}) 
#' of the top and bottom of ranked lists, and the Jaccard index (J) of selected genes.
#' 
#' @author Fabian Schmich
#' @export
#' @seealso \code{\link{Phenotype}}
#' @seealso \code{\link{plot.concordance}}
#' @seealso \code{\link{rbo}}
#' 
#' @param ... The phenotypes to be evaluated for concordance
#' @param min.overlap The minimum number of overlapping genes required
#' @param cor.method A character string indicating which correlation coefficient is to be computed
#' @param cor.use A character string indicating a method for computing correlation in the presence of missing values
#' @param rbo.p The weighting parameter for rank biased overlap (rbo) in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param rbo.k The evaluation depth for rank biased overlap extrapolation
#' @param rbo.mid The mid point to split a ranked list, e.g. in order to split positive and negative scores choose mid=0
#' @return A \code{\link{concordance}} object with the following elements:
#' \item{pair.test}{Indicator of compared phenotypes}
#' \item{cor}{The correlation between pairs of phenotypes for the same gene}
#' \item{rbo.top}{The rank biased overlap of genes evaluated at the top of the ranked list}
#' \item{rbo.bottom}{The rank biased overlap of genes evaluated at the bottom of the ranked list}
#' \item{jaccard}{The Jaccard index of selected genes}
concordance <- function(..., min.overlap=1, cor.method="spearman", cor.use="pairwise.complete.obs", rbo.p=0.95, rbo.k=NULL, rbo.mid=NULL) {
  phenotypes <- list(...)
  if(all(sapply(phenotypes, is, class2="Phenotype")) == FALSE) {
    stop("Not all elements for comparison are Phenotype objects.")
  }
  # Find all pairs of phenotypes
  pairs <- expand.grid(1:length(phenotypes), 1:length(phenotypes))
  pairs <- pairs[which(pairs[,1] < pairs[,2]),]
  cors <- rbo.top <- rbo.bottom <- jaccard <- vector(length=nrow(pairs))
  for (r in 1:nrow(pairs)) {
    p1 <- scores(phenotypes[[pairs[r,1]]])
    p2 <- scores(phenotypes[[pairs[r,2]]])
    if (length(intersect(names(p1), names(p2))) < min.overlap) {
      warn(sprintf("Minimal overlap not met by pair %d - %d.", pairs[r,1], pairs[r,2]))
      cors[r] <- rbo.top[r] <- rbo.bottom[r] <- jaccard[r] <- NA
    } else {
      uids <- unique(names(p1), names(p2))
      p1 <- p1[match(uids, names(p1))]
      p2 <- p2[match(uids, names(p2))]
      cors[r] <- cor(x=p1, y=p2, use=cor.use, method=cor.method)
      rbo.top[r] <- rbo(s=p1, t=p2, side="top", p=rbo.p, mid=rbo.mid)
      rbo.bottom[r] <- rbo(s=p1, t=p2, side="bottom", p=rbo.p, mid=rbo.mid)
      jaccard[r] <- .jaccard(names(p1), names(p2))    
    }
  }
  res <- data.frame(test.pair=paste(pairs[,1], pairs[,2], sep="-"),
                    cor=cors,
                    rbo.top=rbo.top,
                    rbo.bottom=rbo.bottom,
                    jaccard=jaccard)
  class(res) <- "concordance"
  return(res)
}

#' Plot concordance 
#' 
#' Plots boxplots of concordance evaluated between multiple Phenotype objects. Measures include the correlation (rho) between pairs of
#' phenotypes for the same gene, the rank biased overlap (rbo) of the top and bottom of ranked lists, and the
#' Jaccard index (J) of selected genes.
#' 
#' @author Fabian Schmich
#' @export
#' 
#' @param x The data of class \code{\link{concordance}}
plot.concordance <- function(x) {
  if (!require(ggplot2) | !require(reshape2)) {
    warning("You may want to install ggplot2 and reshape2 for prettier plots.")
    boxplot(x[-1])
  } else {
    x <- melt(data.frame(x[1:5]), id.vars=c("test.pair"), variable.name="measure")
    x$measure <- factor(x$measure, levels=c("cor", "rbo.top", "rbo.bottom", "jaccard"))
    ggplot(data=x, aes(x=measure, y=value, colour=measure)) + 
      geom_boxplot(outlier.size=0, width=0.85, na.rm=T) +
      #       scale_x_discrete(labels=c(expression(rho), expression(rbo["" %down% ""]), expression(rbo["" %up% ""]), "J")) +
      xlab("") + ylab("") +
      theme_bw() +
      theme(legend.position="none",
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"))  
  }
}

#' Rank biased overlap (Webber et al., 2010)
#' 
#' Evaluates the rank biased overlap (rbo) of two ranked lists based on formula based on (32) from 
#' "A Similarity Measure for Indefinite Rankings" (Webber et al.). Two ranked lists with high rbo are
#' very similar, wheras low rbo indicates dissimilar lists. rbo ranges between 0 and 1. In this method
#' the extrapolated version of rbo is implemented.
#' 
#' @author Fabian Schmich
#' @export
#' @seealso \code{\link{concordance}}
#' 
#' @param s List 1
#' @param t List 2
#' @param p Weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k Evaluation depth for extrapolation
#' @param side Evaluate similarity between the top or the bottom of the ranked lists
#' @param mid Set the mid point to for example only consider positive or negative scores
#' @return rank biased overlap (rbo)
rbo <- function(s, t, p, k=floor(max(length(s), length(t))/2), side=c("top", "bottom"), mid=NULL) {
  side <- match.arg(side)
  if (!is.numeric(s) | !is.numeric(t))
    stop("Input vectors are not numeric.")
  if (is.null(names(s)) | is.null(names(t)))
    stop("Input vectors are not named.")
  ids <- switch(side,
                "top"=list(s=.select.ids(s, "top", mid), t=.select.ids(t, "top", mid)),
                "bottom"=list(s=.select.ids(s, "bottom", mid), t=.select.ids(t, "bottom", mid))
  )
  min(1, .rbo.ext.unequal.length(ids$s, ids$t, p, k))
}

#' Select top or bottom names of ranked vector
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param x The ranked list
#' @param side The side to be evaluated ("top" or "bottom" of ranked list)
#' @param mid The mid point to split a list, e.g. to split between positive and negative values choose mid=0
#' @return A vector of selected identifiers
.select.ids <- function(x, side=c("top", "bottom"), mid=NULL) {
  side <- match.arg(side)
  if (side == "top")  {
    x <- sort(x, decreasing=TRUE)
    if (is.null(mid))
      return(names(x))
    else 
      return(names(x)[which(x > mid)])
  } else if (side == "bottom") {
    x <- sort(x, decreasing=FALSE)
    if (is.null(mid)) 
      return(names(x))
    else 
      return(names(x)[which(x < mid)])
  }
}

#' Rank biased overlap formula based on (32) from "A Similarity Measure for Indefinite Rankings" (Webber et al.)
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param x List 1
#' @param y List 2
#' @param p The weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k The evaluation depth
#' @return The rank biased overlap between x and y
.rbo.ext.unequal.length <- function(x, y, p, k) {
  if (length(x) <= length(y)) {
    S <- x
    L <- y
  } else {
    S <- y
    L <- x
  }
  l <- min(k, length(L))
  s <- min(k, length(S))
  
  Xd <- sapply(1:l, function(i) length(intersect(S[1:i], L[1:i])))
  ((1-p) / p) *
    ((sum(Xd / seq(1, l) * p^seq(1, l))) +
       (sum(Xd[s] * (seq(s+1, l) - s) / (s * seq(s+1, l)) * p^seq(s+1, l)))) +
    ((Xd[l] - Xd[s]) / l + (Xd[s] / s)) * p^l 
}

#' Jaccard index between two sets
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param x Set 1
#' @param y Set 2
#' @return The Jaccard index
.jaccard <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}
