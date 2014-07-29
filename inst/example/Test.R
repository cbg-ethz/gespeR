# Testing things

x <- TargetRelations("~/scratch/gespeR2/package/gespeR/AMBION_K.rds")
#x <- unloadMatrix(x)
y <- Phenotypes("~/scratch/gespeR2/package/gespeR/AMBION_K.txt", type="SSP", col.id=1, col.score=2)

y <- y[sample(1500)]

f <- join(x, y)

z <- f$targets
#z <- unloadMatrix(z, writeMatrix=FALSE)
z <- unloadMatrix(z, writeMatrix=TRUE)
path(z) <- "~/test.rds"
z <- unloadMatrix(z, writeMatrix=TRUE)
z <- unloadMatrix(z,overwrite=T)
z <- loadMatrix(z)

g.stab <- gespeR(phenotypes=f$phenotypes, target.relations=f$targets, mode="stability", ncores=7, EV=10, threshold=0.6)
g.cv <- gespeR(phenotypes=f$phenotypes, target.relations=f$targets, mode="cv", ncores=5, alpha=0.5)



p <- readRDS("~/scratch/gespeR2/Data/Phenotypes/normCHTS_cn5_ic5/cellHTS2/BRUCELLA_QU-G_B.rds")
y <- Phenotypes(p$scored, channel="cells.all", sample=1)

