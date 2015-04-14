# Code for OU paper
# Gavin Thomas 2014

library(TESS)
library(motmot)
library(caper)

# setwd to location of project

out <- as.data.frame(matrix(NA, ncol=7, nrow=56000))

for (i in 1:1000) {

	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0,MRCA=TRUE)[[1]]
    tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/yule25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}	
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[i,1] <- "yule25"
	out[i,2] <- i
	out[i,3] <- bm_mod$logLikelihood
	out[i,4] <- ou_mod$MaximumLikelihood
	out[i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[i,6] <- ou_mod$Alpha
	out[i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0,MRCA=TRUE)[[1]]
    tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/yule50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)	
	out[1000+i,1] <- "yule50"
	out[1000+i,2] <- i
	out[1000+i,3] <- bm_mod$logLikelihood
	out[1000+i,4] <- ou_mod$MaximumLikelihood
	out[1000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[1000+i,6] <- ou_mod$Alpha
	out[1000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/yule100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[2000+i,1] <- "yule100"
	out[2000+i,2] <- i
	out[2000+i,3] <- bm_mod$logLikelihood
	out[2000+i,4] <- ou_mod$MaximumLikelihood
	out[2000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[2000+i,6] <- ou_mod$Alpha
	out[2000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/yule150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[3000+i,1] <- "yule150"
	out[3000+i,2] <- i
	out[3000+i,3] <- bm_mod$logLikelihood
	out[3000+i,4] <- ou_mod$MaximumLikelihood
	out[3000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[3000+i,6] <- ou_mod$Alpha
	out[3000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/yule200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[4000+i,1] <- "yule200"
	out[4000+i,2] <- i
	out[4000+i,3] <- bm_mod$logLikelihood
	out[4000+i,4] <- ou_mod$MaximumLikelihood
	out[4000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[4000+i,6] <- ou_mod$Alpha
	out[4000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/yule500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[5000+i,1] <- "yule500"
	out[5000+i,2] <- i
	out[5000+i,3] <- bm_mod$logLikelihood
	out[5000+i,4] <- ou_mod$MaximumLikelihood
	out[5000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[5000+i,6] <- ou_mod$Alpha
	out[5000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/yule1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/yule1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[6000+i,1] <- "yule1000"
	out[6000+i,2] <- i
	out[6000+i,3] <- bm_mod$logLikelihood
	out[6000+i,4] <- ou_mod$MaximumLikelihood
	out[6000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[6000+i,6] <- ou_mod$Alpha
	out[6000+i,7] <- mean(ft[[1]][,"I.prime"])

	

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[7000+i,1] <- "bdLow25"
	out[7000+i,2] <- i
	out[7000+i,3] <- bm_mod$logLikelihood
	out[7000+i,4] <- ou_mod$MaximumLikelihood
	out[7000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[7000+i,6] <- ou_mod$Alpha
	out[7000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[8000+i,1] <- "bdLow50"
	out[8000+i,2] <- i
	out[8000+i,3] <- bm_mod$logLikelihood
	out[8000+i,4] <- ou_mod$MaximumLikelihood
	out[8000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[8000+i,6] <- ou_mod$Alpha
	out[8000+i,7] <- mean(ft[[1]][,"I.prime"])

	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[9000+i,1] <- "bdLow100"
	out[9000+i,2] <- i
	out[9000+i,3] <- bm_mod$logLikelihood
	out[9000+i,4] <- ou_mod$MaximumLikelihood
	out[9000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[9000+i,6] <- ou_mod$Alpha
	out[9000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[10000+i,1] <- "bdLow150"
	out[10000+i,2] <- i
	out[10000+i,3] <- bm_mod$logLikelihood
	out[10000+i,4] <- ou_mod$MaximumLikelihood
	out[10000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[10000+i,6] <- ou_mod$Alpha
	out[10000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[11000+i,1] <- "bdLow200"
	out[11000+i,2] <- i
	out[11000+i,3] <- bm_mod$logLikelihood
	out[11000+i,4] <- ou_mod$MaximumLikelihood
	out[11000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[11000+i,6] <- ou_mod$Alpha
	out[11000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[12000+i,1] <- "bdLow500"
	out[12000+i,2] <- i
	out[12000+i,3] <- bm_mod$logLikelihood
	out[12000+i,4] <- ou_mod$MaximumLikelihood
	out[12000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[12000+i,6] <- ou_mod$Alpha
	out[12000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0.25,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdLow1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdLow1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[13000+i,1] <- "bdLow1000"
	out[13000+i,2] <- i
	out[13000+i,3] <- bm_mod$logLikelihood
	out[13000+i,4] <- ou_mod$MaximumLikelihood
	out[13000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[13000+i,6] <- ou_mod$Alpha
	out[13000+i,7] <- mean(ft[[1]][,"I.prime"])

	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[14000+i,1] <- "bdMid25"
	out[14000+i,2] <- i
	out[14000+i,3] <- bm_mod$logLikelihood
	out[14000+i,4] <- ou_mod$MaximumLikelihood
	out[14000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[14000+i,6] <- ou_mod$Alpha
	out[14000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[15000+i,1] <- "bdMid50"
	out[15000+i,2] <- i
	out[15000+i,3] <- bm_mod$logLikelihood
	out[15000+i,4] <- ou_mod$MaximumLikelihood
	out[15000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[15000+i,6] <- ou_mod$Alpha
	out[15000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[16000+i,1] <- "bdMid100"
	out[16000+i,2] <- i
	out[16000+i,3] <- bm_mod$logLikelihood
	out[16000+i,4] <- ou_mod$MaximumLikelihood
	out[16000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[16000+i,6] <- ou_mod$Alpha
	out[16000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[17000+i,1] <- "bdMid150"
	out[17000+i,2] <- i
	out[17000+i,3] <- bm_mod$logLikelihood
	out[17000+i,4] <- ou_mod$MaximumLikelihood
	out[17000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[17000+i,6] <- ou_mod$Alpha
	out[17000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[18000+i,1] <- "bdMid200"
	out[18000+i,2] <- i
	out[18000+i,3] <- bm_mod$logLikelihood
	out[18000+i,4] <- ou_mod$MaximumLikelihood
	out[18000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[18000+i,6] <- ou_mod$Alpha
	out[18000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[19000+i,1] <- "bdMid500"
	out[19000+i,2] <- i
	out[19000+i,3] <- bm_mod$logLikelihood
	out[19000+i,4] <- ou_mod$MaximumLikelihood
	out[19000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[19000+i,6] <- ou_mod$Alpha
	out[19000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0.5,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdMid1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdMid1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[20000+i,1] <- "bdMid1000"
	out[20000+i,2] <- i
	out[20000+i,3] <- bm_mod$logLikelihood
	out[20000+i,4] <- ou_mod$MaximumLikelihood
	out[20000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[20000+i,6] <- ou_mod$Alpha
	out[20000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	
	
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[21000+i,1] <- "bdHigh25"
	out[21000+i,2] <- i
	out[21000+i,3] <- bm_mod$logLikelihood
	out[21000+i,4] <- ou_mod$MaximumLikelihood
	out[21000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[21000+i,6] <- ou_mod$Alpha
	out[21000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[22000+i,1] <- "bdHigh50"
	out[22000+i,2] <- i
	out[22000+i,3] <- bm_mod$logLikelihood
	out[22000+i,4] <- ou_mod$MaximumLikelihood
	out[22000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[22000+i,6] <- ou_mod$Alpha
	out[22000+i,7] <- mean(ft[[1]][,"I.prime"])


	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[23000+i,1] <- "bdHigh100"
	out[23000+i,2] <- i
	out[23000+i,3] <- bm_mod$logLikelihood
	out[23000+i,4] <- ou_mod$MaximumLikelihood
	out[23000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[23000+i,6] <- ou_mod$Alpha
	out[23000+i,7] <- mean(ft[[1]][,"I.prime"])


	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[24000+i,1] <- "bdHigh150"
	out[24000+i,2] <- i
	out[24000+i,3] <- bm_mod$logLikelihood
	out[24000+i,4] <- ou_mod$MaximumLikelihood
	out[24000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[24000+i,6] <- ou_mod$Alpha
	out[24000+i,7] <- mean(ft[[1]][,"I.prime"])


	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[25000+i,1] <- "bdHigh200"
	out[25000+i,2] <- i
	out[25000+i,3] <- bm_mod$logLikelihood
	out[25000+i,4] <- ou_mod$MaximumLikelihood
	out[25000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[25000+i,6] <- ou_mod$Alpha
	out[25000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[26000+i,1] <- "bdHigh500"
	out[26000+i,2] <- i
	out[26000+i,3] <- bm_mod$logLikelihood
	out[26000+i,4] <- ou_mod$MaximumLikelihood
	out[26000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[26000+i,6] <- ou_mod$Alpha
	out[26000+i,7] <- mean(ft[[1]][,"I.prime"])

	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0.75,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/bdHigh1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/bdHigh1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[27000+i,1] <- "bdHigh1000"
	out[27000+i,2] <- i
	out[27000+i,3] <- bm_mod$logLikelihood
	out[27000+i,4] <- ou_mod$MaximumLikelihood
	out[27000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[27000+i,6] <- ou_mod$Alpha
	out[27000+i,7] <- mean(ft[[1]][,"I.prime"])

	
	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[28000+i,1] <- "slowIncrease25"
	out[28000+i,2] <- i
	out[28000+i,3] <- bm_mod$logLikelihood
	out[28000+i,4] <- ou_mod$MaximumLikelihood
	out[28000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[28000+i,6] <- ou_mod$Alpha	
	out[28000+i,7] <- mean(ft[[1]][,"I.prime"])

	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[29000+i,1] <- "slowIncrease50"
	out[29000+i,2] <- i
	out[29000+i,3] <- bm_mod$logLikelihood
	out[29000+i,4] <- ou_mod$MaximumLikelihood
	out[29000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[29000+i,6] <- ou_mod$Alpha
	out[29000+i,7] <- mean(ft[[1]][,"I.prime"])

	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[30000+i,1] <- "slowIncrease100"
	out[30000+i,2] <- i
	out[30000+i,3] <- bm_mod$logLikelihood
	out[30000+i,4] <- ou_mod$MaximumLikelihood
	out[30000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[30000+i,6] <- ou_mod$Alpha
	out[30000+i,7] <- mean(ft[[1]][,"I.prime"])

	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[31000+i,1] <- "slowIncrease150"
	out[31000+i,2] <- i
	out[31000+i,3] <- bm_mod$logLikelihood
	out[31000+i,4] <- ou_mod$MaximumLikelihood
	out[31000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[31000+i,6] <- ou_mod$Alpha
	out[31000+i,7] <- mean(ft[[1]][,"I.prime"])

	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[32000+i,1] <- "slowIncrease200"
	out[32000+i,2] <- i
	out[32000+i,3] <- bm_mod$logLikelihood
	out[32000+i,4] <- ou_mod$MaximumLikelihood
	out[32000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[32000+i,6] <- ou_mod$Alpha
	out[32000+i,7] <- mean(ft[[1]][,"I.prime"])

	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[33000+i,1] <- "slowIncrease500"
	out[33000+i,2] <- i
	out[33000+i,3] <- bm_mod$logLikelihood
	out[33000+i,4] <- ou_mod$MaximumLikelihood
	out[33000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[33000+i,6] <- ou_mod$Alpha
	out[33000+i,7] <- mean(ft[[1]][,"I.prime"])

	l <- function(x) { return(x^2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowIncrease1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowIncrease1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[34000+i,1] <- "slowIncrease1000"
	out[34000+i,2] <- i
	out[34000+i,3] <- bm_mod$logLikelihood
	out[34000+i,4] <- ou_mod$MaximumLikelihood
	out[34000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[34000+i,6] <- ou_mod$Alpha
	out[34000+i,7] <- mean(ft[[1]][,"I.prime"])

	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[35000+i,1] <- "rapidIncrease25"
	out[35000+i,2] <- i
	out[35000+i,3] <- bm_mod$logLikelihood
	out[35000+i,4] <- ou_mod$MaximumLikelihood
	out[35000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[35000+i,6] <- ou_mod$Alpha	
	out[35000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[36000+i,1] <- "rapidIncrease50"
	out[36000+i,2] <- i
	out[36000+i,3] <- bm_mod$logLikelihood
	out[36000+i,4] <- ou_mod$MaximumLikelihood
	out[36000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[36000+i,6] <- ou_mod$Alpha
	out[36000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[37000+i,1] <- "rapidIncrease100"
	out[37000+i,2] <- i
	out[37000+i,3] <- bm_mod$logLikelihood
	out[37000+i,4] <- ou_mod$MaximumLikelihood
	out[37000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[37000+i,6] <- ou_mod$Alpha
	out[37000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[38000+i,1] <- "rapidIncrease150"
	out[38000+i,2] <- i
	out[38000+i,3] <- bm_mod$logLikelihood
	out[38000+i,4] <- ou_mod$MaximumLikelihood
	out[38000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[38000+i,6] <- ou_mod$Alpha
	out[38000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[39000+i,1] <- "rapidIncrease200"
	out[39000+i,2] <- i
	out[39000+i,3] <- bm_mod$logLikelihood
	out[39000+i,4] <- ou_mod$MaximumLikelihood
	out[39000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[39000+i,6] <- ou_mod$Alpha
	out[39000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[40000+i,1] <- "rapidIncrease500"
	out[40000+i,2] <- i
	out[40000+i,3] <- bm_mod$logLikelihood
	out[40000+i,4] <- ou_mod$MaximumLikelihood
	out[40000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[40000+i,6] <- ou_mod$Alpha
	out[40000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidIncrease1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidIncrease1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[41000+i,1] <- "rapidIncrease1000"
	out[41000+i,2] <- i
	out[41000+i,3] <- bm_mod$logLikelihood
	out[41000+i,4] <- ou_mod$MaximumLikelihood
	out[41000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[41000+i,6] <- ou_mod$Alpha
	out[41000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[42000+i,1] <- "slowDecrease25"
	out[42000+i,2] <- i
	out[42000+i,3] <- bm_mod$logLikelihood
	out[42000+i,4] <- ou_mod$MaximumLikelihood
	out[42000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[42000+i,6] <- ou_mod$Alpha	
	out[42000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[43000+i,1] <- "slowDecrease50"
	out[43000+i,2] <- i
	out[43000+i,3] <- bm_mod$logLikelihood
	out[43000+i,4] <- ou_mod$MaximumLikelihood
	out[43000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[43000+i,6] <- ou_mod$Alpha
	out[43000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[44000+i,1] <- "slowDecrease100"
	out[44000+i,2] <- i
	out[44000+i,3] <- bm_mod$logLikelihood
	out[44000+i,4] <- ou_mod$MaximumLikelihood
	out[44000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[44000+i,6] <- ou_mod$Alpha
	out[44000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[45000+i,1] <- "slowDecrease150"
	out[45000+i,2] <- i
	out[45000+i,3] <- bm_mod$logLikelihood
	out[45000+i,4] <- ou_mod$MaximumLikelihood
	out[45000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[45000+i,6] <- ou_mod$Alpha
	out[45000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[46000+i,1] <- "slowDecrease200"
	out[46000+i,2] <- i
	out[46000+i,3] <- bm_mod$logLikelihood
	out[46000+i,4] <- ou_mod$MaximumLikelihood
	out[46000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[46000+i,6] <- ou_mod$Alpha
	out[46000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[47000+i,1] <- "slowDecrease500"
	out[47000+i,2] <- i
	out[47000+i,3] <- bm_mod$logLikelihood
	out[47000+i,4] <- ou_mod$MaximumLikelihood
	out[47000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[47000+i,6] <- ou_mod$Alpha
	out[47000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.5) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/slowDecrease1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/slowDecrease1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[48000+i,1] <- "slowDecrease1000"
	out[48000+i,2] <- i
	out[48000+i,3] <- bm_mod$logLikelihood
	out[48000+i,4] <- ou_mod$MaximumLikelihood
	out[48000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[48000+i,6] <- ou_mod$Alpha
	out[48000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidDecrease25.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease25.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease25.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[49000+i,1] <- "rapidDecrease25"
	out[49000+i,2] <- i
	out[49000+i,3] <- bm_mod$logLikelihood
	out[49000+i,4] <- ou_mod$MaximumLikelihood
	out[49000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[49000+i,6] <- ou_mod$Alpha	
	out[49000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidDecrease50.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease50.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease50.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[50000+i,1] <- "rapidDecrease50"
	out[50000+i,2] <- i
	out[50000+i,3] <- bm_mod$logLikelihood
	out[50000+i,4] <- ou_mod$MaximumLikelihood
	out[50000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[50000+i,6] <- ou_mod$Alpha
	out[50000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidDecrease100.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease100.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease100.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[51000+i,1] <- "rapidDecrease100"
	out[51000+i,2] <- i
	out[51000+i,3] <- bm_mod$logLikelihood
	out[51000+i,4] <- ou_mod$MaximumLikelihood
	out[51000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[51000+i,6] <- ou_mod$Alpha
	out[51000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidDecrease150.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease150.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease150.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[52000+i,1] <- "rapidDecrease150"
	out[52000+i,2] <- i
	out[52000+i,3] <- bm_mod$logLikelihood
	out[52000+i,4] <- ou_mod$MaximumLikelihood
	out[52000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[52000+i,6] <- ou_mod$Alpha
	out[52000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidDecrease200.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease200.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease200.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[53000+i,1] <- "rapidDecrease200"
	out[53000+i,2] <- i
	out[53000+i,3] <- bm_mod$logLikelihood
	out[53000+i,4] <- ou_mod$MaximumLikelihood
	out[53000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[53000+i,6] <- ou_mod$Alpha
	out[53000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
     tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
    write.tree(tr, file="output/rapidDecrease500.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease500.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease500.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[54000+i,1] <- "rapidDecrease500"
	out[54000+i,2] <- i
	out[54000+i,3] <- bm_mod$logLikelihood
	out[54000+i,4] <- ou_mod$MaximumLikelihood
	out[54000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[54000+i,6] <- ou_mod$Alpha
	out[54000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	l <- function(x) { return(x^0.2) }
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,l,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
    tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
	tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/rapidDecrease1000.tre", append=TRUE)
	dat <- transformPhylo.sim(tr, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease1000.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/rapidDecrease1000.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[55000+i,1] <- "rapidDecrease1000"
	out[55000+i,2] <- i
	out[55000+i,3] <- bm_mod$logLikelihood
	out[55000+i,4] <- ou_mod$MaximumLikelihood
	out[55000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[55000+i,6] <- ou_mod$Alpha
	out[55000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	
		
	print(i)
	write.csv(out, file="output/OU_simulations.csv")
	save(out, file="OU_simulations.rda")
}

###### Simulating trees with error
library(TESS)
library(motmot)
library(caper)

setwd("~/Google Drive/CurrentResearch/OrnsteinUhlenbeckDiatribe")

out <- as.data.frame(matrix(NA, ncol=7, nrow=77000))

for (i in 1:1000) {

	
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
    tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	
	write.tree(tr, file="output/error/yule25_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule25_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule25_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE, )
    sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
    out[i,1] <- "yule25_err1"
	out[i,2] <- i
	out[i,3] <- bm_mod$logLikelihood
	out[i,4] <- ou_mod$MaximumLikelihood
	out[i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[i,6] <- ou_mod$Alpha
	out[i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule50_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule50_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule50_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[1000+i,1] <- "yule50_err1"
	out[1000+i,2] <- i
	out[1000+i,3] <- bm_mod$logLikelihood
	out[1000+i,4] <- ou_mod$MaximumLikelihood
	out[1000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[1000+i,6] <- ou_mod$Alpha
	out[1000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule100_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule100_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule100_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[2000+i,1] <- "yule100_err1"
	out[2000+i,2] <- i
	out[2000+i,3] <- bm_mod$logLikelihood
	out[2000+i,4] <- ou_mod$MaximumLikelihood
	out[2000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[2000+i,6] <- ou_mod$Alpha
	out[2000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule150_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule150_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule150_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[3000+i,1] <- "yule150_err1"
	out[3000+i,2] <- i
	out[3000+i,3] <- bm_mod$logLikelihood
	out[3000+i,4] <- ou_mod$MaximumLikelihood
	out[3000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[3000+i,6] <- ou_mod$Alpha
	out[3000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule200_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule200_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule200_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[4000+i,1] <- "yule200_err1"
	out[4000+i,2] <- i
	out[4000+i,3] <- bm_mod$logLikelihood
	out[4000+i,4] <- ou_mod$MaximumLikelihood
	out[4000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[4000+i,6] <- ou_mod$Alpha
	out[4000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule500_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule500_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule500_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[5000+i,1] <- "yule500_err1"
	out[5000+i,2] <- i
	out[5000+i,3] <- bm_mod$logLikelihood
	out[5000+i,4] <- ou_mod$MaximumLikelihood
	out[5000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[5000+i,6] <- ou_mod$Alpha
	out[5000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule1000_err1.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.01*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule1000_err1.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule1000_err1.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[6000+i,1] <- "yule1000_err1"
	out[6000+i,2] <- i
	out[6000+i,3] <- bm_mod$logLikelihood
	out[6000+i,4] <- ou_mod$MaximumLikelihood
	out[6000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[6000+i,6] <- ou_mod$Alpha
	out[6000+i,7] <- mean(ft[[1]][,"I.prime"])

	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule25_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule25_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule25_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE, )
    sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
    out[7000+i,1] <- "yule25_err5"
	out[7000+i,2] <- i
	out[7000+i,3] <- bm_mod$logLikelihood
	out[7000+i,4] <- ou_mod$MaximumLikelihood
	out[7000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[7000+i,6] <- ou_mod$Alpha
	out[7000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule50_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule50_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule50_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[8000+i,1] <- "yule50_err5"
	out[8000+i,2] <- i
	out[8000+i,3] <- bm_mod$logLikelihood
	out[8000+i,4] <- ou_mod$MaximumLikelihood
	out[8000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[8000+i,6] <- ou_mod$Alpha
	out[8000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule100_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule100_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule100_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[9000+i,1] <- "yule100_err5"
	out[9000+i,2] <- i
	out[9000+i,3] <- bm_mod$logLikelihood
	out[9000+i,4] <- ou_mod$MaximumLikelihood
	out[9000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[9000+i,6] <- ou_mod$Alpha
	out[9000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule150_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule150_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule150_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[10000+i,1] <- "yule150_err5"
	out[10000+i,2] <- i
	out[10000+i,3] <- bm_mod$logLikelihood
	out[10000+i,4] <- ou_mod$MaximumLikelihood
	out[10000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[10000+i,6] <- ou_mod$Alpha
	out[10000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule200_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule200_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule200_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[11000+i,1] <- "yule200_err5"
	out[11000+i,2] <- i
	out[11000+i,3] <- bm_mod$logLikelihood
	out[11000+i,4] <- ou_mod$MaximumLikelihood
	out[11000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[11000+i,6] <- ou_mod$Alpha
	out[11000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule500_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule500_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule500_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[12000+i,1] <- "yule500_err5"
	out[12000+i,2] <- i
	out[12000+i,3] <- bm_mod$logLikelihood
	out[12000+i,4] <- ou_mod$MaximumLikelihood
	out[12000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[12000+i,6] <- ou_mod$Alpha
	out[12000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule1000_err5.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.05*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule1000_err5.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule1000_err5.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[13000+i,1] <- "yule1000_err5"
	out[13000+i,2] <- i
	out[13000+i,3] <- bm_mod$logLikelihood
	out[13000+i,4] <- ou_mod$MaximumLikelihood
	out[13000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[13000+i,6] <- ou_mod$Alpha
	out[13000+i,7] <- mean(ft[[1]][,"I.prime"])

	

	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=26,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule25_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule25_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule25_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE, )
    sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
    out[14000+i,1] <- "yule25_err10"
	out[14000+i,2] <- i
	out[14000+i,3] <- bm_mod$logLikelihood
	out[14000+i,4] <- ou_mod$MaximumLikelihood
	out[14000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[14000+i,6] <- ou_mod$Alpha
	out[14000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=51,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule50_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule50_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule50_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[15000+i,1] <- "yule50_err10"
	out[15000+i,2] <- i
	out[15000+i,3] <- bm_mod$logLikelihood
	out[15000+i,4] <- ou_mod$MaximumLikelihood
	out[15000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[15000+i,6] <- ou_mod$Alpha
	out[15000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=101,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule100_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule100_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule100_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[16000+i,1] <- "yule100_err10"
	out[16000+i,2] <- i
	out[16000+i,3] <- bm_mod$logLikelihood
	out[16000+i,4] <- ou_mod$MaximumLikelihood
	out[16000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[16000+i,6] <- ou_mod$Alpha
	out[16000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=151,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule150_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule150_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule150_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[17000+i,1] <- "yule150_err10"
	out[17000+i,2] <- i
	out[17000+i,3] <- bm_mod$logLikelihood
	out[17000+i,4] <- ou_mod$MaximumLikelihood
	out[17000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[17000+i,6] <- ou_mod$Alpha
	out[17000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=201,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule200_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule200_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule200_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[18000+i,1] <- "yule200_err10"
	out[18000+i,2] <- i
	out[18000+i,3] <- bm_mod$logLikelihood
	out[18000+i,4] <- ou_mod$MaximumLikelihood
	out[18000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[18000+i,6] <- ou_mod$Alpha
	out[18000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=501,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule500_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule500_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule500_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[19000+i,1] <- "yule500_err10"
	out[19000+i,2] <- i
	out[19000+i,3] <- bm_mod$logLikelihood
	out[19000+i,4] <- ou_mod$MaximumLikelihood
	out[19000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[19000+i,6] <- ou_mod$Alpha
	out[19000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	tr <- sim.globalBiDe.taxa(n=1,nTaxa=1001,max=10,1,0,MRCA=TRUE)[[1]]
	tipidx <- which(tr$edge.length==min(tr$edge.length[match(1:Ntip(tr), tr$edge[,2])]))
	tip.d <- tr$tip.label[tr$edge[tipidx[1],2]]
	tr <- drop.tip(tr, tip.d)
    tr$tip.label <- paste("t", 1:Ntip(tr), sep="")
	write.tree(tr, file="output/error/yule1000_err10.tre", append=TRUE)
	tr2 <- tr
	tr2$edge.length[match(1:Ntip(tr), tr$edge[,2])] <- tr$edge.length[match(1:Ntip(tr), tr$edge[,2])] + 0.1*max(node.depth.edgelength(tr))
	dat <- transformPhylo.sim(tr2, n=1, model="bm")
	if (i==1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule1000_err10.txt", append=FALSE, col.names=TRUE, row.names=FALSE)}
	if (i!=1) {write.table(t(dat[order(rownames(dat)),]), file="output/error/yule1000_err10.txt", append=TRUE, col.names=FALSE, row.names=FALSE)}
	bm_mod <- transformPhylo.ML(dat, tr, model="bm", modelCIs=FALSE)
	ou_mod <- transformPhylo.ML(dat, tr, model="OU", modelCIs=FALSE)	
	sprich <- data.frame(species=tr$tip.label, rich=rep(1, Ntip(tr)))
	cdat <- comparative.data(phy=tr, dat=sprich, names.col=species)
	ft <- fusco.test(cdat, rich=rich, reps=1)
	out[20000+i,1] <- "yule1000_err10"
	out[20000+i,2] <- i
	out[20000+i,3] <- bm_mod$logLikelihood
	out[20000+i,4] <- ou_mod$MaximumLikelihood
	out[20000+i,5] <- 2*(ou_mod$MaximumLikelihood - bm_mod$logLikelihood)
	out[20000+i,6] <- ou_mod$Alpha
	out[20000+i,7] <- mean(ft[[1]][,"I.prime"])
	
	print(i)
	write.csv(out, file="output/error/OU_error_simulations.csv")
	save(out, file="output/error/OU_error_simulations.rda")
	}








load("OU_simulations.rda")
out <- out2
out <- data.frame(type=out[,1], index=as.numeric(out[,2]), bm=as.numeric(out[,3]), ou=as.numeric(out[,4]), LR=as.numeric(out[,5]), alpha=as.numeric(out[,6]), I.prime=as.numeric(out[,7]))


t1err <- function(x) { return(x > 3.84) }
HPD <- function(x) {qu <- quantile(x, probs=seq(0,1,0.025))
								   return(c(qu[2],qu[40]))
								   }


xx1 <- grep("yule", out[,1])
xx2 <- grep("err", out[,1])
xx <- out[setdiff(xx1, xx2), ]
yule <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yule_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleerr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("bdLow", out[,1]), ]
bdlow <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdlow_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdlowerr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("bdMid", out[,1]), ]
bdmid <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdmid_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdmiderr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("bdHigh", out[,1]), ]
bdhigh <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdhigh_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdhigherr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("slowIncrease", out[,1]), ]
slowIncrease <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
slowIncrease_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
slowIncreaseerr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("rapidIncrease", out[,1]), ]
rapidIncrease <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
rapidIncrease_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
rapidIncreaseerr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("slowDecrease", out[,1]), ]
slowDecrease <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
slowDecrease_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
slowDecreaseerr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("rapidDecrease", out[,1]), ]
rapidDecrease <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
rapidDecrease_HPD <- aggregate(xx[,6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
rapidDecreaseerr <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx1 <- grep("_err1", out[,1])
xx2 <- grep("_err10", out[,1])
xx <- out[setdiff(xx1, xx2), ]
yuleErr1 <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr1_sd <- aggregate(xx[,5:6], by=list(xx[,1]), sd)[c(5,6,1,3,4,7,2),]
yuleErr1_err <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("_err5", out[,1]), ]
yuleErr5 <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr5_sd <- aggregate(xx[,5:6], by=list(xx[,1]), sd)[c(5,6,1,3,4,7,2),]
yuleErr5_err <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("_err10", out[,1]), ]
yuleErr10 <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr10_sd <- aggregate(xx[,5:6], by=list(xx[,1]), sd)[c(5,6,1,3,4,7,2),]
yuleErr10_err <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

### Get MCMC results
mcmcres <- read.csv("output/MCMC_OU01Uni.csv")

xx1 <- grep("yule", mcmcres[,1])
xx2 <- grep("err", mcmcres[,1])
xx <- mcmcres[setdiff(xx1, xx2), ]
yulemcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
yulemcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
yulemcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]


xx <- mcmcres[grep("bdLow", mcmcres[,1]),]
bdlowmcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
bdlowmcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
bdlowmcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]

xx <- mcmcres[grep("bdMid", mcmcres[,1]),]
bdmidmcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
bdmidmcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
bdmidmcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]

xx <- mcmcres[grep("bdHigh", mcmcres[,1]),]
bdhighmcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
bdhighmcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
bdhighmcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]

xx <- mcmcres[grep("slowIncrease", mcmcres[,1]),]
slowIncreasemcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
slowIncreasemcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
slowIncreasemcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]

xx <- mcmcres[grep("rapidIncrease", mcmcres[,1]),]
rapidIncreasemcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
rapidIncreasemcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
rapidIncreasemcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]

xx <- mcmcres[grep("slowDecrease", mcmcres[,1]),]
slowDecreasemcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
slowDecreasemcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
slowDecreasemcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]

xx <- mcmcres[grep("rapidDecrease", mcmcres[,1]),]
rapidDecreasemcmc <- aggregate(xx[,"OU"], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
rapidDecreasemcmcLHPD <- aggregate(xx[,"X2.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]
rapidDecreasemcmcUHPD <- aggregate(xx[,"X97.50."], by=list(xx[,"tree"]), mean)[c(5,6,1,3,4,7,2),]




### Tree imbalance ###

pdf(file="figs/OU_imbalace.pdf", width=6, height=6)
plot(out$alpha ~ out$I.prime, pch=16, cex=0.5, col=rgb(0,0,0,0.5), ylab="alpha", xlab="Imbalance (I)")
abline(v=0.5, col="red")
dev.off()


### Alpha is fucked plots ###

colscheme <- list(rgb(166,97,26, maxColorValue=255), rgb(223,194,125, maxColorValue=255), rgb(128,205,193, maxColorValue=255), rgb(1,133,113, maxColorValue=255)) 

pdf(file="figs/OU_tree_sizeshape_mean.pdf", width=8, height=12)
	par(mfrow=c(3,2), mai=c(0.8,0.8,0.5,0.1))
	
	# Constant rates
	plot(1:7, yuleerr/1000, ylim=c(0,0.12), pch=16, xlab="Tree size", ylab="Type I error", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]])
	points(1:7, bdlowerr/1000, pch=16, col=colscheme[[2]])
	points(1:7, bdmiderr/1000, pch=16, col=colscheme[[3]])
	points(1:7, bdhigherr/1000, pch=16, col=colscheme[[4]])
	
	mtext("(A) Constant rate trees", 3, adj=0, line=0.2) 
	legend(x=5, y=0.12, c("d/b=0", "d/b=0.25","d/b=0.5","d/b=0.75"), pch=16, bty="n", col=unlist(colscheme))
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.05, col="grey", lty="dashed")
	
	plot(0.9:6.9, yule[,3], ylim=c(0,6.5), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]])
	points(0.9:6.9, as.matrix(yule_HPD[,2])[,1], pch="-", col=colscheme[[1]])
	points(0.9:6.9, as.matrix(yule_HPD[,2])[,2], pch="-", col=colscheme[[1]])
	points(0.9:6.9, bdlow[,3], pch=16, col=colscheme[[2]])
	points(0.9:6.9, as.matrix(bdlow_HPD[,2])[,1], pch="-", col=colscheme[[2]])
	points(0.9:6.9, as.matrix(bdlow_HPD[,2])[,2], pch="-", col=colscheme[[2]])
	points(0.9:6.9, bdmid[,3], pch=16, col=colscheme[[3]])
	points(0.9:6.9, as.matrix(bdmid_HPD[,2])[,1], pch="-", col=colscheme[[3]])
	points(0.9:6.9, as.matrix(bdmid_HPD[,2])[,2], pch="-", col=colscheme[[3]])
	points(0.9:6.9, bdhigh[,3], pch=16, col=colscheme[[4]])
	points(0.9:6.9, as.matrix(bdhigh_HPD[,2])[,1], pch="-", col=colscheme[[4]])
	points(0.9:6.9, as.matrix(bdhigh_HPD[,2])[,2], pch="-", col=colscheme[[4]])

points(1.1:7.1, yulemcmc[,2], pch=15, col=colscheme[[1]])
points(1.1:7.1, yulemcmcLHPD[,2], pch="-", col=colscheme[[1]])
points(1.1:7.1, yulemcmcUHPD[,2], pch="-", col=colscheme[[1]])

points(1.1:7.1, bdlowmcmc[,2], pch=15, col=colscheme[[2]])
points(1.1:7.1, bdlowmcmcLHPD[,2], pch="-", col=colscheme[[2]])
points(1.1:7.1, bdlowmcmcUHPD[,2], pch="-", col=colscheme[[2]])

points(1.1:7.1, bdmidmcmc[,2], pch=15, col=colscheme[[3]])
points(1.1:7.1, bdmidmcmcLHPD[,2], pch="-", col=colscheme[[3]])
points(1.1:7.1, bdmidmcmcUHPD[,2], pch="-", col=colscheme[[3]])

points(1.1:7.1, bdhighmcmc[,2], pch=15, col=colscheme[[4]])
points(1.1:7.1, bdhighmcmcLHPD[,2], pch="-", col=colscheme[[4]])
points(1.1:7.1, bdhighmcmcUHPD[,2], pch="-", col=colscheme[[4]])

	mtext("(B) Constant rate trees", 3, adj=0, line=0.2)
	legend(x=5, y=6.5, c("d/b=0", "d/b=0.25","d/b=0.5","d/b=0.75"), pch=16, bty="n", col=unlist(colscheme))
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.0, col="grey", lty="dashed")	
	
	# Increasing rates
	plot(1:7, slowIncreaseerr/1000, ylim=c(0,0.2), pch=16, xlab="Tree size", ylab="Type I error", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]])
	points(1:7, rapidIncreaseerr/1000, pch=16, col=colscheme[[3]])

	mtext("(C) Tippy trees", 3, adj=0, line=0.2) 
	legend(x=5, y=0.2, c("Slow increase", "Fast increase"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]))
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.05, col="grey", lty="dashed")

	plot(0.9:6.9, slowIncrease[,3], ylim=c(0,7), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]])
	points(0.9:6.9, as.matrix(slowIncrease_HPD[,2])[,1], pch="-", col=colscheme[[1]])
	points(0.9:6.9, as.matrix(slowIncrease_HPD[,2])[,2], pch="-", col=colscheme[[1]])
	points(0.9:6.9, rapidIncrease[,3], pch=16, col=colscheme[[3]])
	points(0.9:6.9, as.matrix(rapidIncrease_HPD[,2])[,1], pch="-", col=colscheme[[3]])
	points(0.9:6.9, as.matrix(rapidIncrease_HPD[,2])[,2], pch="-", col=colscheme[[3]])


points(1.1:7.1, slowIncreasemcmc[,2], pch=15, col=colscheme[[1]])
points(1.1:7.1, slowIncreasemcmcLHPD[,2], pch="-", col=colscheme[[1]])
points(1.1:7.1, slowIncreasemcmcUHPD[,2], pch="-", col=colscheme[[1]])

points(1.1:7.1, rapidIncreasemcmc[,2], pch=15, col=colscheme[[3]])
points(1.1:7.1, rapidIncreasemcmcLHPD[,2], pch="-", col=colscheme[[3]])
points(1.1:7.1, rapidIncreasemcmcUHPD[,2], pch="-", col=colscheme[[3]])


	mtext("(D) Tippy trees", 3, adj=0, line=0.2) 
	legend(x=5, y=7, c("Slow increase", "Fast increase"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]))
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.0, col="grey", lty="dashed")
	
	
	# Decreasing rates
	plot(1:7, slowDecreaseerr/1000, ylim=c(0,0.12), pch=16, xlab="Tree size", ylab="Type I error", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]])
	points(1:7, rapidDecreaseerr/1000, pch=16, col=colscheme[[3]])

	mtext("(E) Rooty trees", 3, adj=0, line=0.2) 
	legend(x=5, y=0.12, c("Slow decrease", "Fast decrease"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]))
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.05, col="grey", lty="dashed")
	
	plot(0.9:6.9, slowDecrease[,3], ylim=c(0,6.5), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5,  bty="l", col=colscheme[[1]])
	points(0.9:6.9, as.matrix(slowDecrease_HPD[,2])[,1], pch="-", col=colscheme[[1]])
	points(0.9:6.9, as.matrix(slowDecrease_HPD[,2])[,2], pch="-", col=colscheme[[1]])
	points(0.9:6.9, rapidDecrease[,3], pch=16, tcl=0.5, col=colscheme[[3]])
	points(0.9:6.9, as.matrix(rapidDecrease_HPD[,2])[,1], pch="-", col=colscheme[[3]])
	points(0.9:6.9, as.matrix(rapidDecrease_HPD[,2])[,2], pch="-", col=colscheme[[3]])

points(1.1:7.1, slowDecreasemcmc[,2], pch=15, col=colscheme[[1]])
points(1.1:7.1, slowDecreasemcmcLHPD[,2], pch="-", col=colscheme[[1]])
points(1.1:7.1, slowDecreasemcmcUHPD[,2], pch="-", col=colscheme[[1]])

points(1.1:7.1, rapidDecreasemcmc[,2], pch=15, col=colscheme[[3]])
points(1.1:7.1, rapidDecreasemcmcLHPD[,2], pch="-", col=colscheme[[3]])
points(1.1:7.1, rapidDecreasemcmcUHPD[,2], pch="-", col=colscheme[[3]])


	mtext("(F) Rooty trees", 3, adj=0, line=0.2) 
	legend(x=5, y=6.5, c("Slow decrease", "Fast decrease"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]))
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.0, col="grey", lty="dashed")
	
dev.off()



pdf(file="figs/OU_measurement_error.pdf", width=8, height=4)
	par(mfrow=c(1,2), mai=c(0.8,0.8,0.5,0.1))
	
	plot(1:7, yuleErr1_err/1000, ylim=c(0,1), pch=1, xlab="Tree size", ylab="Proportion favouring OU", xaxt="n", tcl=0.5, bty="l")
	points(1:7, yuleErr5_err/1000, pch=15)
	points(1:7, yuleErr10_err/1000, pch=16)
	
	mtext("(A)", 3, adj=0, line=0.2) 
	legend(x=5, y=0.4, c("1% error", "5% error", "10% error"), pch=c(1,15,16), bty="n")
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.05, col="red", lty="dashed")

	plot(1:7, yuleErr1[,3], ylim=c(0,4), pch=1, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l")
	points(1:7, yuleErr5[,3], pch=15)
	points(1:7, yuleErr10[,3], pch=16)
	
	mtext("(B)", 3, adj=0, line=0.2) 
	legend(x=5, y=4.2, c("1% error", "5% error", "10% error"), pch=c(1,15,16), bty="n")	
	axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
	abline(h=0.0, col="red", lty="dashed")

dev.off()



pdf(file="output/OUsims_typeI.pdf", width=12, height=8)
par(mfrow=c(3,4))
plot(1:7, yuleerr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Yule", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, bdlowerr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Birth-death, ext frac=0.25", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, bdmiderr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Birth-death, ext frac=0.5", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, bdhigherr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Birth-death, ext frac=0.75", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, slowIncreaseerr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Increasing rates (slow)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, rapidIncreaseerr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Increasing rates (fast)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, slowDecreaseerr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Decreasing rates (slow)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, rapidDecreaseerr/1000, ylim=c(0,0.16), pch=16, xlab="Tree size", ylab="Type I error", main="Decreasing rates (fast)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, yuleErr1_err/1000, ylim=c(0,1), pch=16, xlab="Tree size", ylab="Type I error", main="Yule 1% error", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, yuleErr5_err/1000, ylim=c(0,1), pch=16, xlab="Tree size", ylab="Type I error", main="Yule 5% error", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, yuleErr10_err/1000, ylim=c(0,1), pch=16, xlab="Tree size", ylab="Type I error", main="Yule 10% error", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

dev.off()



pdf(file="output/OUsims_alpha.pdf", width=12, height=8)
par(mfrow=c(3,4))
plot(1:7, yule[,3], ylim=c(0,0.3), pch=16, xlab="Tree size", ylab="alpha", main="Yule", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, bdlow[,3], ylim=c(0,0.3), pch=16, xlab="Tree size", ylab="alpha", main="Birth-death, ext frac=0.25", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, bdmid[,3], ylim=c(0,0.3), pch=16, xlab="Tree size", ylab="alpha", main="Birth-death, ext frac=0.5", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, bdhigh[,3], ylim=c(0,0.3), pch=16, xlab="Tree size", ylab="alpha", main="Birth-death, ext frac=0.75", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, slowIncrease[,3], ylim=c(0,1.5), pch=16, xlab="Tree size", ylab="alpha", main="Increasing rates (slow)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, rapidIncrease[,3], ylim=c(0,1.5), pch=16, xlab="Tree size", ylab="alpha", main="Increasing rates (fast)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, slowDecrease[,3], ylim=c(0,1.5), pch=16, xlab="Tree size", ylab="alpha", main="Decreasing rates (slow)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, rapidDecrease[,3], ylim=c(0,1.5), pch=16, xlab="Tree size", ylab="alpha", main="Decreasing rates (fast)", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.0, col="grey")

plot(1:7, yuleErr1[,3], ylim=c(0,4), pch=16, xlab="Tree size", ylab="alpha", main="Yule 1% error", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, yuleErr5[,3], ylim=c(0,4), pch=16, xlab="Tree size", ylab="alpha", main="Yule 5% error", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

plot(1:7, yuleErr10[,3], ylim=c(0,4), pch=16, xlab="Tree size", ylab="alpha", main="Yule 10% error", xaxt="n")
axis(side=1, at=1:7, labels=c(25,50,100,150,200,500,1000))
abline(h=0.05, col="grey")

dev.off()







###### Likelihood profiles ########
# Using trees of 100 tips.        #
# Three trees.                    #
# One rapid increase.             #
# One rapid decrease.             #
# One yule.                       #
###################################

setwd("~/Google Drive/CurrentResearch/OrnsteinUhlenbeckDiatribe")

trsTippy <- read.tree("output/rapidIncrease50.tre")
trsRooty <- read.tree("output/rapidDecrease50.tre")
trsYule <- read.tree("output/yule50.tre")


datTippy <- read.table("output/rapidIncrease50.txt", header=TRUE)
datRooty <- read.table("output/rapidDecrease50.txt", header=TRUE)
datYule <- read.table("output/yule50.txt", header=TRUE)


# Randomly choose data set
set.seed(120101)
samp <- sample(1:1000, 10, replace=FALSE)


trTippy <- trsTippy[[samp[1]]]
datTippySamp <- t(datTippy[samp[1],])

trRooty <- trsRooty[[samp[2]]]
datRootySamp <- t(datRooty[samp[2],])

trYule <- trsYule[[samp[3]]]
datYuleSamp <- t(datYule[samp[3],])


alphaSamp <- seq(0.0001, 10, 0.01)
lamSamp <- seq(0, 1, 0.001)

outTippyOU <- matrix(NA, ncol=11, nrow=length(alphaSamp))
outTippyOU[,1] <- alphaSamp

outTippyLam <- matrix(NA, ncol=11, nrow=length(lamSamp))
outTippyLam[,1] <- lamSamp


for (j in 1:10) {
	for (i in 1:length(alphaSamp)) {
	
	outTippyOU[i,j+1]  <- transformPhylo.ll(t(datTippy[samp[j],]), trsTippy[[samp[j]]], model="OU", alpha=alphaSamp[i])[[2]]
	outTippyLam[i,j+1]  <- transformPhylo.ll(t(datTippy[samp[j],]), trsTippy[[samp[j]]], model="lambda", lambda= lamSamp[i])[[2]]
	print(i)}
}





outRootyOU <- matrix(NA, ncol=11, nrow=length(alphaSamp))
outRootyOU[,1] <- alphaSamp

outRootyLam <- matrix(NA, ncol=11, nrow=length(lamSamp))
outRootyLam[,1] <- lamSamp


for (j in 1:10) {
	for (i in 1:length(alphaSamp)) {
	
	outRootyOU[i,j+1]  <- transformPhylo.ll(t(datRooty[samp[j],]), trsRooty[[samp[j]]], model="OU", alpha=alphaSamp[i])[[2]]
	outRootyLam[i,j+1]  <- transformPhylo.ll(t(datRooty[samp[j],]), trsRooty[[samp[j]]], model="lambda", lambda= lamSamp[i])[[2]]
	print(i)}
}



outYuleOU <- matrix(NA, ncol=11, nrow=length(alphaSamp))
outYuleOU[,1] <- alphaSamp

outYuleLam <- matrix(NA, ncol=11, nrow=length(lamSamp))
outYuleLam[,1] <- lamSamp


for (j in 1:10) {
	for (i in 1:length(alphaSamp)) {
	
	outYuleOU[i,j+1]  <- transformPhylo.ll(t(datYule[samp[j],]), trsYule[[samp[j]]], model="OU", alpha=alphaSamp[i])[[2]]
	outYuleLam[i,j+1]  <- transformPhylo.ll(t(datYule[samp[j],]), trsYule[[samp[j]]], model="lambda", lambda= lamSamp[i])[[2]]
	print(i)}
}


write.csv(outTippyOU, file="output/ProfileDataTippyOU50.csv", row.names=FALSE)
write.csv(outTippyLam, file="output/ProfileDataTippyLam50.csv", row.names=FALSE)
write.csv(outRootyOU, file="output/ProfileDataRootyOU50.csv", row.names=FALSE)
write.csv(outRootyLam, file="output/ProfileDataRootyLam50.csv", row.names=FALSE)
write.csv(outYuleOU, file="output/ProfileDataYuleOU50.csv", row.names=FALSE)
write.csv(outYuleLam, file="output/ProfileDataYuleLam50.csv", row.names=FALSE)







outRooty <- rep(NA, length(alphaSamp))
outYule <- rep(NA, length(alphaSamp))






outRooty[i]  <- transformPhylo.ll(datRootySamp, trRooty, model="OU", alpha=alphaSamp[i])[[2]]
outYule[i]  <- transformPhylo.ll(datYuleSamp, trYule, model="OU", alpha=alphaSamp[i])[[2]]
print(i)}

pdf("figs/ProfilePlots.pdf", height=30, width=10)
par(mfcol=c(10,3))
plot(outTippy ~alphaSamp, type="l", col="red", ylab="log-likelihood", xlab="alpha", main="Tippy tree")
plot(outRooty ~alphaSamp, type="l", col="red", ylab="log-likelihood", xlab="alpha", main="Rooty tree")
plot(outYule ~alphaSamp, type="l", col="red", ylab="log-likelihood", xlab="alpha", main="Yule tree")
dev.off()


write.csv(cbind(alpha=alphaSamp, tippyLnL=outTippy, rootyLnL=outRooty, yuleLnL=outYule), file="output/ProfileData.csv")


