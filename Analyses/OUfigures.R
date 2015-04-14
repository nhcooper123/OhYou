###################################################
#### Plot proportion error and alpha estimates ####
###################################################
setwd("/Users/ghthomas/Google Drive/CurrentResearch/OrnsteinUhlenbeckDiatribe")


load("OU_simulations.rda")
out <- data.frame(type=out[,1], index=as.numeric(out[,2]), bm=as.numeric(out[,3]), ou=as.numeric(out[,4]), LR=as.numeric(out[,5]), alpha=as.numeric(out[,6]), I.prime=as.numeric(out[,7]))


t1err <- function(x) { return(x > 3.84) }
HPD <- function(x) {qu <- quantile(x, probs=seq(0,1,0.025))
    return(c(qu[2],qu[40]))
}

BF <- function (x) { return(x > 2) }


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



### Get MCMC results
mcmcres <- read.csv("data/MCMC_ExpLam10.csv")

xx1 <- grep("yule", mcmcres[,1])
xx2 <- grep("err", mcmcres[,1])
xx <- mcmcres[setdiff(xx1, xx2), ]
yulemcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yulemcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yulemcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("bdLow", mcmcres[,1]),]
bdlowmcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdlowmcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdlowmcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("bdMid", mcmcres[,1]),]
bdmidmcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdmidmcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdmidmcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("bdHigh", mcmcres[,1]),]
bdhighmcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdhighmcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdhighmcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("slowIncrease", mcmcres[,1]),]
slowIncreasemcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
slowIncreasemcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
slowIncreasemcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("rapidIncrease", mcmcres[,1]),]
rapidIncreasemcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
rapidIncreasemcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
rapidIncreasemcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("slowDecrease", mcmcres[,1]),]
slowDecreasemcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
slowDecreasemcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
slowDecreasemcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("rapidDecrease", mcmcres[,1]),]
rapidDecreasemcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
rapidDecreasemcmcLHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
rapidDecreasemcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])






colscheme <- list(rgb(166,97,26, maxColorValue=255), rgb(223,194,125, maxColorValue=255), rgb(128,205,193, maxColorValue=255), rgb(1,133,113, maxColorValue=255))

pdf(file="figs/OU_tree_sizeshape_mean_ExpLam10_27Feb.pdf", width=12, height=16)
par(mfrow=c(3,2), mai=c(0.4,0.4,0.4,0.1), mgp=c(1.5,0.25,0))

# Constant rates
plotseq <- seq(1,9,by=1.25)
plot(plotseq-0.4, yuleerr/1000, ylim=c(0,0.12), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
points(plotseq-0.3, bdlowerr/1000, pch=16, col=colscheme[[2]], cex=1.25)
points(plotseq-0.2, bdmiderr/1000, pch=16, col=colscheme[[3]], cex=1.25)
points(plotseq-0.1, bdhigherr/1000, pch=16, col=colscheme[[4]], cex=1.25)

points(plotseq+0.1, yulemcmc_err/1000, ylim=c(0,0.12), pch=18, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25)
points(plotseq+0.2, bdlowmcmc_err/1000, pch=18, col=colscheme[[2]], cex=1.25)
points(plotseq+0.3, bdmidmcmc_err/1000, pch=18, col=colscheme[[3]], cex=1.25)
points(plotseq+0.4, bdhighmcmc_err/1000, pch=18, col=colscheme[[4]], cex=1.25)


mtext("(A) Constant speciation rate", 3, adj=0, line=0.2)
legend(x=7, y=0.12, c("d/b=0", "d/b=0.25","d/b=0.5","d/b=0.75"), pch=16, bty="n", col=unlist(colscheme), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.4, yule[,3], ylim=c(0,1.5), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
for (i in 1:length(yule[,3])) {
    lines(x=c((plotseq[i]-0.4),(plotseq[i]-0.4)), y=c(yule[i,3], as.matrix(yule_HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.4),(plotseq[i]-0.4)), y=c(yule[i,3], as.matrix(yule_HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.3, bdlow[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[2]], cex=1.25)
for (i in 1:length(bdlow[,3])) {
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(bdlow[i,3], as.matrix(bdlow_HPD[,2])[i,1]), col=colscheme[[2]], lty=3)
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(bdlow[i,3], as.matrix(bdlow_HPD[,2])[i,2]), col=colscheme[[2]], lty=3)
}

points(plotseq-0.2, bdmid[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(bdmid[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(bdmid[i,3], as.matrix(bdmid_HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(bdmid[i,3], as.matrix(bdmid_HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}

points(plotseq-0.1, bdhigh[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[4]], cex=1.25)
for (i in 1:length(bdhigh[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(bdhigh[i,3], as.matrix(bdhigh_HPD[,2])[i,1]), col=colscheme[[4]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(bdhigh[i,3], as.matrix(bdhigh_HPD[,2])[i,2]), col=colscheme[[4]], lty=3)
}

#####

points(plotseq+0.1, yulemcmc[,3], ylim=c(0,14), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25)
for (i in 1:length(yulemcmc[,3])) {
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(yulemcmc[i,3], as.matrix(yulemcmcLHPD[,2])[i,1]), col=colscheme[[1]], lty=2)
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(yulemcmc[i,3], as.matrix(yulemcmcLHPD[,2])[i,2]), col=colscheme[[1]], lty=2)
}

points(plotseq+0.2, bdlowmcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[2]], cex=1.25)
for (i in 1:length(bdlowmcmc[,3])) {
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(bdlowmcmc[i,3], as.matrix(bdlowmcmcLHPD[,2])[i,1]), col=colscheme[[2]], lty=2)
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(bdlowmcmc[i,3], as.matrix(bdlowmcmcLHPD[,2])[i,2]), col=colscheme[[2]], lty=2)
}

points(plotseq+0.3, bdmidmcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(bdmidmcmc[,3])) {
    lines(x=c((plotseq[i]+0.3),(plotseq[i]+0.3)), y=c(bdmidmcmc[i,3], as.matrix(bdmidmcmcLHPD[,2])[i,1]), col=colscheme[[3]], lty=2)
    lines(x=c((plotseq[i]+0.3),(plotseq[i]+0.3)), y=c(bdmidmcmc[i,3], as.matrix(bdmidmcmcLHPD[,2])[i,2]), col=colscheme[[3]], lty=2)
}

points(plotseq+0.4, bdhighmcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[4]], cex=1.25)
for (i in 1:length(bdhighmcmc[,3])) {
    lines(x=c((plotseq[i]+0.4),(plotseq[i]+0.4)), y=c(bdhighmcmc[i,3], as.matrix(bdhighmcmcLHPD[,2])[i,1]), col=colscheme[[4]], lty=2)
    lines(x=c((plotseq[i]+0.4),(plotseq[i]+0.4)), y=c(bdhighmcmc[i,3], as.matrix(bdhighmcmcLHPD[,2])[i,2]), col=colscheme[[4]], lty=2)
}


mtext("(B) Constant speciation rate", 3, adj=0, line=0.2)
plotseq <- seq(1,9,by=1.25)
legend(x=7, y=1.5, c("d/b=0", "d/b=0.25","d/b=0.5","d/b=0.75"), pch=16, bty="n", col=unlist(colscheme), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.0, col="grey", lty="dashed")


##########################
#### Increasing rates ####
plot(plotseq-0.2, slowIncreaseerr/1000, ylim=c(0,0.2), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
points(plotseq-0.1, rapidIncreaseerr/1000, pch=16, col=colscheme[[3]], cex=1.25)

points(plotseq+0.1, slowIncreasemcmc_err/1000, ylim=c(0,0.12), pch=18, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25)
points(plotseq+0.2, rapidIncreasemcmc_err/1000, pch=18, col=colscheme[[3]], cex=1.25)


mtext("(C) Accelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=0.2, c("Slow increase", "Fast increase"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.2, slowIncrease[,3], ylim=c(0,7), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
for (i in 1:length(slowIncrease[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowIncrease[i,3], as.matrix(slowIncrease_HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowIncrease[i,3], as.matrix(slowIncrease_HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.1, rapidIncrease[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(rapidIncrease[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidIncrease[i,3], as.matrix(rapidIncrease_HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidIncrease[i,3], as.matrix(rapidIncrease_HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}


points(plotseq+0.1, slowIncreasemcmc[,3], ylim=c(0,1.5), xlim=c(0.5,9), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25)
for (i in 1:length(slowIncreasemcmc[,3])) {
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(slowIncreasemcmc[i,3], as.matrix(slowIncreasemcmcLHPD[,2])[i,1]), col=colscheme[[1]], lty=2)
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(slowIncreasemcmc[i,3], as.matrix(slowIncreasemcmcLHPD[,2])[i,2]), col=colscheme[[1]], lty=2)
}

points(plotseq+0.2, rapidIncreasemcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(rapidIncreasemcmc[,3])) {
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(rapidIncreasemcmc[i,3], as.matrix(rapidIncreasemcmcLHPD[,2])[i,1]), col=colscheme[[3]], lty=2)
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(rapidIncreasemcmc[i,3], as.matrix(rapidIncreasemcmcLHPD[,2])[i,2]), col=colscheme[[3]], lty=2)
}

mtext("(D) Accelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=7, c("Slow increase", "Fast increase"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.0, col="grey", lty="dashed")


##########################
#### Decreasing rates ####
plot(plotseq-0.2, slowDecreaseerr/1000, ylim=c(0,0.12), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
points(plotseq-0.1, rapidDecreaseerr/1000, pch=16, col=colscheme[[3]], cex=1.25)

points(plotseq+0.1, slowDecreasemcmc_err/1000, ylim=c(0,0.12), pch=18, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25)
points(plotseq+0.2, rapidDecreasemcmc_err/1000, pch=18, col=colscheme[[3]], cex=1.25)


mtext("(E) Decelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=0.12, c("Slow decrease", "Fast decrease"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.2, slowDecrease[,3], ylim=c(0,2), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
for (i in 1:length(slowDecrease[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowDecrease[i,3], as.matrix(slowDecrease_HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowDecrease[i,3], as.matrix(slowDecrease_HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.1, rapidDecrease[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(rapidDecrease[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidDecrease[i,3], as.matrix(rapidDecrease_HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidDecrease[i,3], as.matrix(rapidDecrease_HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}


points(plotseq+0.1, slowDecreasemcmc[,3], ylim=c(0,1.5), xlim=c(0.5,9), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25)
for (i in 1:length(slowDecreasemcmc[,3])) {
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(slowDecreasemcmc[i,3], as.matrix(slowDecreasemcmcLHPD[,2])[i,1]), col=colscheme[[1]], lty=2)
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(slowDecreasemcmc[i,3], as.matrix(slowDecreasemcmcLHPD[,2])[i,2]), col=colscheme[[1]], lty=2)
}

points(plotseq+0.2, rapidDecreasemcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(rapidDecreasemcmc[,3])) {
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(rapidDecreasemcmc[i,3], as.matrix(rapidDecreasemcmcLHPD[,2])[i,1]), col=colscheme[[3]], lty=2)
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(rapidDecreasemcmc[i,3], as.matrix(rapidDecreasemcmcLHPD[,2])[i,2]), col=colscheme[[3]], lty=2)
}


mtext("(F) Decelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=2, c("Slow decrease", "Fast decrease"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.0, col="grey", lty="dashed")



dev.off()






###################################################
#### Plots for data with simulated error       ####
###################################################

setwd("/Users/ghthomas/Google Drive/CurrentResearch/OrnsteinUhlenbeckDiatribe")


load("output/error/OU_error_simulations.rda")
out <- data.frame(type=out[,1], index=as.numeric(out[,2]), bm=as.numeric(out[,3]), ou=as.numeric(out[,4]), LR=as.numeric(out[,5]), alpha=as.numeric(out[,6]), I.prime=as.numeric(out[,7]))


t1err <- function(x) { return(x > 3.84) }
HPD <- function(x) {qu <- quantile(x, probs=seq(0,1,0.025))
    return(c(qu[2],qu[40]))
}

BF <- function (x) { return(x > 2) }



xx1 <- grep("_err1", out[,1])
xx2 <- grep("_err10", out[,1])
xx <- out[setdiff(xx1, xx2), ]
yuleErr1 <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr1HPD <- aggregate(xx[,5:6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr1_err <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("_err5", out[,1]), ]
yuleErr5 <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr5HPD <- aggregate(xx[,5:6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr5_err <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("_err10", out[,1]), ]
yuleErr10 <- aggregate(xx[,5:6], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr10HPD <- aggregate(xx[,5:6], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr10_err <- rowSums(aggregate(xx[,5], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])



mcmcres <- read.csv("data/MCMC_ExpLam10.csv")

xx1 <- grep("_err1", mcmcres[,1])
xx2 <- grep("_err10", mcmcres[,1])
xx <- mcmcres[setdiff(xx1, xx2), ]
yuleErr1mcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr1mcmcHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr1mcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("_err5", mcmcres[,1]),]
yuleErr5mcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr5mcmcHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr5mcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])

xx <- mcmcres[grep("_err10", mcmcres[,1]),]
yuleErr10mcmc <- aggregate(xx[,c("Exp_L10_2Bfactor", "Exp_L10_Alpha_Mode")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr10mcmcHPD <- aggregate(xx[,"Exp_L10_Alpha_Mode"], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr10mcmc_err <- rowSums(aggregate(xx[,"Exp_L10_2Bfactor"], by=list(xx[,1]), BF)[c(5,6,1,3,4,7,2),2])






colscheme <- list(rgb(166,97,26, maxColorValue=255), rgb(223,194,125, maxColorValue=255), rgb(128,205,193, maxColorValue=255), rgb(1,133,113, maxColorValue=255))

pdf(file="figs/OUError_tree_sizeshape_mean_Lam10_27Feb.pdf", width=12, height=16/3)
par(mfrow=c(1,2), mai=c(0.5,0.5,0.4,0.1), mgp=c(1.5,0.25,0))

# Constant rates
plotseq <- seq(1,9,by=1.25)
plot(plotseq-0.4, yuleErr1_err/1000, ylim=c(0,1), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1)
points(plotseq-0.3, yuleErr5_err/1000, pch=16, col=colscheme[[2]], cex=1)
points(plotseq-0.2, yuleErr10_err/1000, pch=16, col=colscheme[[3]], cex=1)

points(plotseq+0.1, yuleErr1mcmc_err/1000, ylim=c(0,0.12), pch=18, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1)
points(plotseq+0.2, yuleErr5mcmc_err/1000, pch=18, col=colscheme[[2]], cex=1)
points(plotseq+0.3, yuleErr10mcmc_err/1000, pch=18, col=colscheme[[3]], cex=1)


mtext("(A) Constant rate + error", 3, adj=0, line=0.2)
legend(x=7.5, y=0.3, c("1% error", "5% error", "10% error"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[2]], colscheme[[3]]))
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.3, yuleErr1[,3], ylim=c(0,20), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1)
for (i in 1:length(yuleErr1[,3])) {
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(yuleErr1[i,3], as.matrix(yuleErr1HPD[,3])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(yuleErr1[i,3], as.matrix(yuleErr1HPD[,3])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.2, yuleErr5[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[2]], cex=1)
for (i in 1:length(yuleErr5[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(yuleErr5[i,3], as.matrix(yuleErr5HPD[,3])[i,1]), col=colscheme[[2]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(yuleErr5[i,3], as.matrix(yuleErr5HPD[,3])[i,2]), col=colscheme[[2]], lty=3)
}

points(plotseq-0.1, yuleErr10[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1)
for (i in 1:length(yuleErr10[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(yuleErr10[i,3], as.matrix(yuleErr10HPD[,3])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(yuleErr10[i,3], as.matrix(yuleErr10HPD[,3])[i,2]), col=colscheme[[3]], lty=3)
}



#####

points(plotseq+0.1, yuleErr1mcmc[,3], ylim=c(0,14), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1)
for (i in 1:length(yuleErr1mcmc[,3])) {
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(yuleErr1mcmc[i,3], as.matrix(yuleErr1mcmcHPD[,2])[i,1]), col=colscheme[[1]], lty=2)
    lines(x=c((plotseq[i]+0.1),(plotseq[i]+0.1)), y=c(yuleErr1mcmc[i,3], as.matrix(yuleErr1mcmcHPD[,2])[i,2]), col=colscheme[[1]], lty=2)
}

points(plotseq+0.2, yuleErr5mcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[2]], cex=1)
for (i in 1:length(yuleErr5mcmc[,3])) {
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(yuleErr5mcmc[i,3], as.matrix(yuleErr5mcmcHPD[,2])[i,1]), col=colscheme[[2]], lty=2)
    lines(x=c((plotseq[i]+0.2),(plotseq[i]+0.2)), y=c(yuleErr5mcmc[i,3], as.matrix(yuleErr5mcmcHPD[,2])[i,2]), col=colscheme[[2]], lty=2)
}

points(plotseq+0.3, yuleErr10mcmc[,3], ylim=c(0,8), pch=18, xlab="Tree size", ylab="alpha", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1)
for (i in 1:length(yuleErr10mcmc[,3])) {
    lines(x=c((plotseq[i]+0.3),(plotseq[i]+0.3)), y=c(yuleErr10mcmc[i,3], as.matrix(yuleErr10mcmcHPD[,2])[i,1]), col=colscheme[[3]], lty=2)
    lines(x=c((plotseq[i]+0.3),(plotseq[i]+0.3)), y=c(yuleErr10mcmc[i,3], as.matrix(yuleErr10mcmcHPD[,2])[i,2]), col=colscheme[[3]], lty=2)
}




mtext("(B) Constant rate + error", 3, adj=0, line=0.2)
plotseq <- seq(1,9,by=1.25)
legend(x=7.5, y=20, c("1% error", "5% error", "10% error"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[2]], colscheme[[3]]))
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
abline(h=0.0, col="grey", lty="dashed")



dev.off()








###################################
### Tree imbalance              ###
###################################
pdf(file="figs/OU_imbalace.pdf", width=6, height=6)
plot(out$alpha ~ out$I.prime, pch=16, cex=0.5, col=rgb(0,0,0,0.5), ylab="alpha", xlab="Imbalance (I)")
abline(v=0.5, col="red")
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


