# Simulation of gespeR example data

require(Matrix)

# Set statics
N <- 1000
p <- 1500
n.effects <- n.effects <- 0.05 * N # number of effect genes
on.strength <- 0.75


# Set beta ground truth
beta.truth <- vector(length=p)
beta.truth[sample(1:p, n.effects, replace=F)] <- rnorm(n.effects, 0, 3)


# Set on-targets
# Each siRNA has one on-target
# Each gene is targeted at least once
on.targets <- vector(length=N)
stopifnot(p >= N)
for (i in 1:N) {
  on.targets[i] <- sample(setdiff(1:p, on.targets), 1)
}

simulate.data <- function(beta.truth, on.targets, N, p, n.effects, on.strength) {
  # Initialise target relations
  X <- Matrix(0, nrow=N, ncol=p, dimnames=list(rownames=sprintf("siRNAID_%.4d", 1:N), colnames=sprintf("geneID_%.4d", 1:p)))
  
  X[cbind(1:N, on.targets)] <- on.strength
  
  # Set off-targets
  # Number of off-targets is Normal(3e-2 * N, 3e-3 * N) distributed
  # Strength of off-target is Beta(2, 5) distributed
  for (i in 1:N) {
    on <- which(X[i,] == on.strength)
    n.off <- floor(rnorm(n=1, mean=3e-2 * N, sd=3e-3 * N))
    X[i,sample(1:p, n.off, replace=F)] <- rbeta(n.off, shape1=2, shape2=5)
    X[i, on] <- on.strength
  }
  # format(object.size(X), units="Mb")
  
  # Set phenotypes
  # Represented are Z-scored phenotypic readouts
  Y <- X %*% beta.truth + rnorm(N, 0, 0.5)
  
  return(list(X=X, Y=Y))
}


# Simulate 4 screens
dat <- lapply(LETTERS[1:4], function(x) {
  simulate.data(beta.truth, on.targets, N, p, n.effects, on.strength)
})
names(dat) <- sprintf("screen_%s", LETTERS[1:4])
names(beta.truth) <- colnames(dat[[1]]$X)
dat <- c(list(beta=beta.truth), dat)




lapply(names(Y), function(x) {
  write.table(data.frame(siRNAID=rownames(Y[[x]]), Scores=Y[[x]][,1]), file=sprintf("~/scratch/gespeR2/Bioconductor/gespeR/inst/extdata/Phenotypes_%s.txt", x), row.names=F, quote=F, sep="\t")
})

lapply(names(X), function(x) {
  saveRDS(X[[x]], file=sprintf("~/scratch/gespeR2/Bioconductor/gespeR/inst/extdata/TR_%s.rds", x))
})

