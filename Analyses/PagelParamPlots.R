setwd("/Users/ghthomas/Google Drive/CurrentResearch/OrnsteinUhlenbeckDiatribe")


out <- read.csv("data/DeltaLambdaKappaMLRes2_Working.csv")

#### Delta
xx1 <- grep("yule", out[,1])
xx2 <- grep("err", out[,1])
xx <- out[setdiff(xx1, xx2), ]
yule <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yule_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleerr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("bdLow", out[,1]), ]
bdlow <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdlow_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdlowerr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("bdMid", out[,1]), ]
bdmid <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdmid_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdmiderr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("bdHigh", out[,1]), ]
bdhigh <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
bdhigh_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
bdhigherr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("slowIncrease", out[,1]), ]
slowIncrease <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
slowIncrease_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
slowIncreaseerr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("rapidIncrease", out[,1]), ]
rapidIncrease <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
rapidIncrease_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
rapidIncreaseerr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("slowDecrease", out[,1]), ]
slowDecrease <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
slowDecrease_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
slowDecreaseerr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("rapidDecrease", out[,1]), ]
rapidDecrease <- aggregate(xx[,c("Delta_LR", "Delta")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
rapidDecrease_HPD <- aggregate(xx[,c("Delta")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
rapidDecreaseerr <- rowSums(aggregate(xx[,"Delta_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])





colscheme <- list(rgb(166,97,26, maxColorValue=255), rgb(223,194,125, maxColorValue=255), rgb(128,205,193, maxColorValue=255), rgb(1,133,113, maxColorValue=255))

pdf(file="figs/Delta_tree_sizeshape_mean__27Feb.pdf", width=12, height=16)
par(mfrow=c(3,2), mai=c(0.4,0.4,0.4,0.1), mgp=c(1.5,0.25,0))

# Constant rates
plotseq <- seq(1,9,by=1.25)
plot(plotseq-0.4, yuleerr/1000, ylim=c(0,0.2), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
points(plotseq-0.3, bdlowerr/1000, pch=16, col=colscheme[[2]], cex=1.25)
points(plotseq-0.2, bdmiderr/1000, pch=16, col=colscheme[[3]], cex=1.25)
points(plotseq-0.1, bdhigherr/1000, pch=16, col=colscheme[[4]], cex=1.25)


mtext("(A) Constant speciation rate", 3, adj=0, line=0.2)
legend(x=7, y=0.20, c("d/b=0", "d/b=0.25","d/b=0.5","d/b=0.75"), pch=16, bty="n", col=unlist(colscheme), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.4, yule[,3], ylim=c(0,4), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
for (i in 1:length(yule[,3])) {
    lines(x=c((plotseq[i]-0.4),(plotseq[i]-0.4)), y=c(yule[i,3], as.matrix(yule_HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.4),(plotseq[i]-0.4)), y=c(yule[i,3], as.matrix(yule_HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.3, bdlow[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[2]], cex=1.25)
for (i in 1:length(bdlow[,3])) {
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(bdlow[i,3], as.matrix(bdlow_HPD[,2])[i,1]), col=colscheme[[2]], lty=3)
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(bdlow[i,3], as.matrix(bdlow_HPD[,2])[i,2]), col=colscheme[[2]], lty=3)
}

points(plotseq-0.2, bdmid[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(bdmid[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(bdmid[i,3], as.matrix(bdmid_HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(bdmid[i,3], as.matrix(bdmid_HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}

points(plotseq-0.1, bdhigh[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[4]], cex=1.25)
for (i in 1:length(bdhigh[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(bdhigh[i,3], as.matrix(bdhigh_HPD[,2])[i,1]), col=colscheme[[4]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(bdhigh[i,3], as.matrix(bdhigh_HPD[,2])[i,2]), col=colscheme[[4]], lty=3)
}


mtext("(B) Constant speciation rate", 3, adj=0, line=0.2)
plotseq <- seq(1,9,by=1.25)
legend(x=7, y=4, c("d/b=0", "d/b=0.25","d/b=0.5","d/b=0.75"), pch=16, bty="n", col=unlist(colscheme), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=1.0, col="grey", lty="dashed")


##########################
#### Increasing rates ####
plot(plotseq-0.2, slowIncreaseerr/1000, ylim=c(0,0.2), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
points(plotseq-0.1, rapidIncreaseerr/1000, pch=16, col=colscheme[[3]], cex=1.25)


mtext("(C) Accelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=0.2, c("Slow increase", "Fast increase"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.2, slowIncrease[,3], ylim=c(0,4), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
for (i in 1:length(slowIncrease[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowIncrease[i,3], as.matrix(slowIncrease_HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowIncrease[i,3], as.matrix(slowIncrease_HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.1, rapidIncrease[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(rapidIncrease[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidIncrease[i,3], as.matrix(rapidIncrease_HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidIncrease[i,3], as.matrix(rapidIncrease_HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}


mtext("(D) Accelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=4, c("Slow increase", "Fast increase"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=1.0, col="grey", lty="dashed")


##########################
#### Decreasing rates ####
plot(plotseq-0.2, slowDecreaseerr/1000, ylim=c(0,0.2), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
points(plotseq-0.1, rapidDecreaseerr/1000, pch=16, col=colscheme[[3]], cex=1.25)


mtext("(E) Decelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=0.20, c("Slow decrease", "Fast decrease"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.2, slowDecrease[,3], ylim=c(0,4), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1.25, cex.axis=1.25, cex.lab=1.25)
for (i in 1:length(slowDecrease[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowDecrease[i,3], as.matrix(slowDecrease_HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(slowDecrease[i,3], as.matrix(slowDecrease_HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.1, rapidDecrease[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Delta", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1.25)
for (i in 1:length(rapidDecrease[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidDecrease[i,3], as.matrix(rapidDecrease_HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(rapidDecrease[i,3], as.matrix(rapidDecrease_HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}



mtext("(F) Decelerating speciation", 3, adj=0, line=0.2)
legend(x=7, y=4, c("Slow decrease", "Fast decrease"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[3]]), cex=1.5)
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5, cex.axis=1.25)
abline(h=1.0, col="grey", lty="dashed")



dev.off()






####### Error plots


xx1 <- grep("_err1", mcmcres[,1])
xx2 <- grep("_err10", mcmcres[,1])
xx <- out[setdiff(xx1, xx2), ]
yuleErr1 <- aggregate(xx[,c("Kappa_LR", "Kappa")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr1HPD <- aggregate(xx[,c("Kappa")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr1_err <- rowSums(aggregate(xx[,"Kappa_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("_err5", mcmcres[,1]),]
yuleErr5 <- aggregate(xx[,c("Kappa_LR", "Kappa")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr5HPD <- aggregate(xx[,c("Kappa")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr5_err <- rowSums(aggregate(xx[,"Kappa_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])

xx <- out[grep("_err10", mcmcres[,1]),]
yuleErr10 <- aggregate(xx[,c("Kappa_LR", "Kappa")], by=list(xx[,1]), mean)[c(5,6,1,3,4,7,2),]
yuleErr10HPD <- aggregate(xx[,c("Kappa")], by=list(xx[,1]), HPD)[c(5,6,1,3,4,7,2),]
yuleErr10_err <- rowSums(aggregate(xx[,"Kappa_LR"], by=list(xx[,1]), t1err)[c(5,6,1,3,4,7,2),2])










colscheme <- list(rgb(166,97,26, maxColorValue=255), rgb(223,194,125, maxColorValue=255), rgb(128,205,193, maxColorValue=255), rgb(1,133,113, maxColorValue=255))

pdf(file="figs/KappaError_tree_sizeshape_mean_27Feb.pdf", width=12, height=16/3)
par(mfrow=c(1,2), mai=c(0.5,0.5,0.4,0.1), mgp=c(1.5,0.25,0))

# Constant rates
plotseq <- seq(1,9,by=1.25)
plot(plotseq-0.4, yuleErr1_err/1000, ylim=c(0,1), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Propn. BM rejected", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1)
points(plotseq-0.3, yuleErr5_err/1000, pch=16, col=colscheme[[2]], cex=1)
points(plotseq-0.2, yuleErr10_err/1000, pch=16, col=colscheme[[3]], cex=1)



mtext("(A) Constant rate + error", 3, adj=0, line=0.2)
legend(x=7.5, y=0.3, c("1% error", "5% error", "10% error"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[2]], colscheme[[3]]))
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
abline(h=0.05, col="grey", lty="dashed")

plot(plotseq-0.3, yuleErr1[,3], ylim=c(0,2), xlim=c(0.5,9), pch=16, xlab="Tree size", ylab="Kappa", xaxt="n", tcl=0.5, bty="l", col=colscheme[[1]], cex=1)
for (i in 1:length(yuleErr1[,3])) {
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(yuleErr1[i,3], as.matrix(yuleErr1HPD[,2])[i,1]), col=colscheme[[1]], lty=3)
    lines(x=c((plotseq[i]-0.3),(plotseq[i]-0.3)), y=c(yuleErr1[i,3], as.matrix(yuleErr1HPD[,2])[i,2]), col=colscheme[[1]], lty=3)
}

points(plotseq-0.2, yuleErr5[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Kappa", xaxt="n", tcl=0.5, bty="l", col=colscheme[[2]], cex=1)
for (i in 1:length(yuleErr5[,3])) {
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(yuleErr5[i,3], as.matrix(yuleErr5HPD[,2])[i,1]), col=colscheme[[2]], lty=3)
    lines(x=c((plotseq[i]-0.2),(plotseq[i]-0.2)), y=c(yuleErr5[i,3], as.matrix(yuleErr5HPD[,2])[i,2]), col=colscheme[[2]], lty=3)
}

points(plotseq-0.1, yuleErr10[,3], ylim=c(0,8), pch=16, xlab="Tree size", ylab="Kappa", xaxt="n", tcl=0.5, bty="l", col=colscheme[[3]], cex=1)
for (i in 1:length(yuleErr10[,3])) {
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(yuleErr10[i,3], as.matrix(yuleErr10HPD[,2])[i,1]), col=colscheme[[3]], lty=3)
    lines(x=c((plotseq[i]-0.1),(plotseq[i]-0.1)), y=c(yuleErr10[i,3], as.matrix(yuleErr10HPD[,2])[i,2]), col=colscheme[[3]], lty=3)
}



mtext("(B) Constant rate + error", 3, adj=0, line=0.2)
plotseq <- seq(1,9,by=1.25)
legend(x=7.5, y=2, c("1% error", "5% error", "10% error"), pch=16, bty="n", col=c(colscheme[[1]], colscheme[[2]], colscheme[[3]]))
axis(side=1, at=plotseq, labels=c(25,50,100,150,200,500,1000), tcl=0.5)
abline(h=1.0, col="grey", lty="dashed")



dev.off()








