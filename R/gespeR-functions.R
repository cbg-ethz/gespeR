#' Bind multiple phenotypes to matrix
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param phenotypes list of phenotypes
#' @return binding of phenotypes (matrix)
.bindphens <- function(phenotypes) {
  if (length(unique(sapply(phenotypes, length))) == 1) {
    do.call("cBind", phenotypes)
  } else {
      warning("Cannot bind phenotypes. Unequal length.")
      return(phenotypes)
  }
}

#' gespeR cross validation
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param SSP The siRNA-specific phenotypes. Single vector for univariate, list of
#' vectors for multivariate phenotypes.
#' @param targets The siRNA-to-gene target relations
#' @param alpha The \code{\link{glmnet}} mixing parameter
#' @param ncores The number of cores for parallel computation
#' @return A list containing the fitted model and used paramers
.gespeR.cv <- function(SSP, targets, alpha, ncores = 1) {
  multivar <- ifelse(ncol(SSP) > 1, TRUE, FALSE)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  model <- cv.glmnet(x = targets,
                     y = as.matrix(SSP),
                     family = ifelse(multivar, "mgaussian", "gaussian"),
                     alpha = alpha,
                     type.measure = "mse",
                     standardize = FALSE,
                     intercept = FALSE,
                     keep = TRUE,
                     parallel = ifelse(ncores > 1, TRUE, FALSE))  
  stopCluster(cl)
  if (multivar) coefficients <- do.call("cBind", coef(model, s = "lambda.1se")) else 
    coefficients <- coef(model, s = "lambda.1se")
  out <- list(type = c("cv"),
              multivar = multivar,
              fit = model$glmnet.fit,
              coefficients = coefficients,
              cv = list(name = model$name,
                      nzero = model$nzero, 
                      cvm = model$cvm, 
                      cvsd = model$cvsd, 
                      cvup = model$cvup, 
                      cvlo = model$cvlo, 
                      foldid = model$foldid, 
                      alpha = alpha),
              stability = list()
  )
  return(out)
}

#' gespeR stability selection
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param SSP The siRNA-specific phenotypes
#' @param targets The siRNA-to-gene target relations
#' @param nboostrap The number of bootstrap samples
#' @param fraction The fraction for each bootstrap sample
#' @param threshold The selection threshold
#' @param EV The expected value of wrongly selected elements
#' @param weakness The weakness parameter for randomised lasso
#' @param ncores The number of cores for parallel computation
#' @return A list containing the fitted model and used paramers
.gespeR.stability <- function(SSP, 
                              targets,
                              nbootstrap = 100,
                              fraction = 0.5, 
                              threshold = 0.75, 
                              EV = 1, 
                              weakness = 1,
                              ncores = 1) {
  stab.out <- stability.selection(x = targets, 
                                  y = SSP,
                                  fraction = fraction, 
                                  threshold = threshold, 
                                  EV = EV, 
                                  nbootstrap = nbootstrap,
                                  weakness = weakness,
                                  ncores = ncores,
                                  intercept = FALSE,
                                  standardize = FALSE)
  out <- list(type = c("stability"),
              multivar = FALSE,
#              fit=stab.out$model, # saving space
              coefficients = stab.out$model$coefficients,
              cv = list(),
              stability = c(stab.out[c("frequency", "selection")], EV = EV, threshold = threshold, q = q, nbootstrap = nbootstrap)
  )
  return(out)
}

#'  Randomized Lasso
#'  
#'  Based on Meinshausen and Buehlmann (2009)
#'  
#'  @author Fabian Schmich
#'  @export
#'  
#'  @param x The design matrix
#'  @param y The response vector
#'  @param weakness The weakness parameter
#'  @param subsample The data subsample (default: none)
#'  @param dfmax The maxiumum number of degrees of freedom
#'  @param lambda The regularisation parameter
#'  @param standardize Indicator, wheter to standardize the design matrix
#'  @param intercept Indicator, whether to fit an intercept
#'  @param ... Additional arguments to \code{\link{glmnet}}
#'  @return A \code{\link{glmnet}} object
#'  @examples
#'  y <- rnorm(50)
#'  x <- matrix(runif(50 * 20), ncol = 20)
#'  lasso.rand(x = x, y = y)
lasso.rand <- function(x, 
                       y,
                       weakness = 1,
                       subsample = 1:nrow(x),
                       dfmax = (ncol(x)+1),
                       lambda = NULL,
                       standardize = FALSE,
                       intercept = FALSE,
                       ...) {
  if (is.matrix(y) && ncol(y) > 1) {
    stop("Randomised lasso not yet implemented for multivariate phenotypes")
  }
  if (is.null(dim(y))) y <- cbind(y)  
  glmnet(x[subsample,], 
         y[subsample,],
         family = "gaussian",
         lambda = lambda,
         penalty.factor = (1 / runif(ncol(x), weakness, 1)),
         alpha = 1,
         dfmax = dfmax,
         standardize = FALSE,
         ...)
}


#' Stability Selection
#' 
#' Based on Meinshausen and Buehlmann (2009)
#' 
#' @author Fabian Schmich
#' 
#' @import doParallel
#' @import foreach
#' @import parallel
#' 
#' @export
#' 
#' @param x The design matrix
#' @param y The response vector
#' @param intercept Indicator, whether to fit an intercept
#' @param nbootstrap The number of bootstrap samples
#' @param fraction The fraction for each bootstrap sample
#' @param threshold The selection threshold
#' @param EV The expected value of wrongly selected elements
#' @param weakness The weakness parameter for randomised lasso
#' @param ncores The number of cores for parallel computation
#' @param ... Additional arguments to \code{\link{lasso.rand}}
#' @return A list containing selected covariates with frequencies, and the fitted model
stability.selection <- function(x, y, 
                                fraction = 0.5, 
                                threshold = 0.75, 
                                EV = 1, 
                                nbootstrap = 100,
                                weakness = 1,
                                intercept = FALSE,
                                ncores = 1,
                                ...) {
  # Dimensions
  n <- nrow(x)
  p <- ncol(x)
  
  # Subsample size
  n.sel <- floor(fraction * n)
  
  # Variable bound
  q <- ceiling(sqrt(EV * p * (2 * threshold - 1))) # (9)
  #   cat(sprintf("\nq = %d\n", q))
  
  # Subsampling
#   registerDoMC(cores=ncores)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  sel.mat <- foreach (b = 1:nbootstrap, .combine =rbind) %dopar% {
    # Current sub-sampled data
    sel <- sample(1:n, n.sel, replace=FALSE)        
    # Get selected model
    fit <- lasso.rand(x=x, y=y, subsample=sel, dfmax=q, weakness=weakness, intercept=intercept, ...)
    return(.select.model(fit, q))
  }
  
  # Get selection frequencies
  freq <- colMeans(sel.mat)
  names(freq) <- colnames(x)  
  sel.current <- which(freq >= threshold)
  names(sel.current) <- colnames(x)[sel.current]
  
  # fit model
  x.sel <- as.matrix(x[,sel.current])
  rownames(x.sel) <- rownames(x)
  colnames(x.sel) <- colnames(x)[sel.current]
  model <- lm(as.matrix(y) ~ as.matrix(x.sel) - 1) # as.matrix() ?
  
  # Stop the cluster
  stopCluster(cl)

  out <- list(model = model, matrix = sel.mat, frequency = freq, selection = sel.current)
  return(out)  
}


#' Model selector (taken from hdi package)
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param fit A (randomised lasso) model
#' @param q The variable bound
#' @return A vector indicating the selected model
.select.model <- function(fit, q) {
  p <- fit$dim[1]
  p.sel <- vector(length=p)
  # Determine non-zero coeficients
  nz <- predict(fit, type="nonzero")
  # Determine largest model that is <= q
  delta <- q - unlist(lapply(nz, length)) # deviation from desired model size
  delta[delta < 0] <- Inf # overshooting not allowed
  nz <- nz[[which.min(delta)]] # takes first occurrence
  p.sel[nz] <- TRUE # set selected to TRUE
  return(p.sel) 
}