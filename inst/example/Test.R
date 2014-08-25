# Testing things

# Load data
# load(system.file("extdata", "simdata.RData", package="gespeR"))
phenos <- sapply(LETTERS[1:4], function(x) {
  Phenotypes(system.file("extdata", sprintf("Phenotypes_screen_%s.txt", x), package="gespeR"), type="SSP", col.id=1, col.score=2)
})

targets <- sapply(LETTERS[1:4], function(x) {
  TargetRelations(system.file("extdata", sprintf("TR_screen_%s.rds", x), package="gespeR"))
})

# Basic functions
print(phenos$C)
plot(phenos$C)
print(targets$C)

# Subsetting and join()
sub.libC.phenos <- phenos$C[sample(500, replace=F)]
sub.libC <- join(targets$C, sub.libC.phenos)

# Unloading and loading values
z <- sub.libC$targets
path(z) <- "/tmp/test.rds"
z <- unloadValues(z, overwrite=TRUE)
print(z)
format(object.size(z), "Kb")
z <- loadValues(z)
format(object.size(z), "Kb")

# Run gespeR in stability mode
ges.run <- lapply(LETTERS[1:4], function(x) {
  gespeR(phenos[[x]], targets[[x]], mode="stability", ncores=2, EV=1, threshold=0.75)
})

# Basic function for gespeR objects
summary(ges.run[[1]])
print(ges.run[[1]])

# Concordance
GSP <- sapply(ges.run, gsp)
conc.gsp <- concordance(GSP[[1]], GSP[[2]], GSP[[3]], GSP[[4]])
conc.ssp <- concordance(phenos[[1]], phenos[[2]], phenos[[3]], phenos[[4]])

plot(conc.ssp)
plot(conc.gsp)