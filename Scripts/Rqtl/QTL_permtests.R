library(qtl)
susc <- read.cross('csv', file='QTL_SNP_map.csv', genotypes=c('AA', 'AB', 'BB'))
summary(susc)

dat <- calc.genoprob(susc, step=1, error.prob=0.01)

set.seed(85842518)
operm2A <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2A, file="perm2A.RData")

set.seed(85842519)
operm2B <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2B, file="perm2B.RData")

set.seed(85842520)
operm2C <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2C, file="perm2C.RData")

set.seed(85842521)
operm2D <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2D, file="perm2D.RData")

set.seed(85842522)
operm2E <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2E, file="perm2E.RData")

set.seed(85842523)
operm2F <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2F, file="perm2F.RData")

set.seed(85842524)
operm2G <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2G, file="perm2G.RData")

set.seed(85842525)
operm2H <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2H, file="perm2H.RData")

set.seed(85842526)
operm2I <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2I, file="perm2I.RData")

set.seed(85842527)
operm2J <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2J, file="perm2J.RData")